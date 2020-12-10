//
// Created by sayan on 6/18/20.
//

#include "hiper3ga.h"
#include "hidx.h"
#include "hfile.h"
#include "metadata.h"
#include "kthread.h"

typedef struct mm_tbuf_s {
    u128_v v_mm; /// query minimizers
    u128_v coef; /// Hough transform coefficient
    u128_v intv; /// intervals on sorted coef
    u32_v reg2mini;

    /// the following are for computing LIS
    uint32_t n, m;
    uint64_t *a;
    size_t *b, *p;

    /// final output
    kvec_t(mm_reg1_t) reg;
    ma_hit_v hits;

    /// bookkeeping
    size_t n_aln;
} mm_tbuf_t;
typedef struct {
    int w, k;
    readinfo_t *readinfo;
    uint64_t block;
    char *buffer;
    void *idx;
    mm_tbuf_t *tbuf;
} shared_data_t;
static inline void mm_tbuf_init(mm_tbuf_t *self) {
    kv_init(self->v_mm); kv_init(self->coef); kv_init(self->intv);
    kv_init(self->reg2mini); kv_init(self->reg);
    self->a = self->b = self->p = 0;
    self->n = self->m = 0;
    kv_init(self->hits);
    self->n_aln = 0;
}
static inline void mm_tbuf_clear(mm_tbuf_t *self) {
    kv_destroy(self->v_mm); kv_destroy(self->coef); kv_destroy(self->intv);
    kv_destroy(self->reg2mini); kv_destroy(self->reg);
    hfree(self->a); hfree(self->b); hfree(self->p);
    kv_destroy(self->hits);
}
static void smp_aligner_fn(void *s, long i, int tid);
static void process_block(shared_data_t *sd, runconf_t *runconf, size_t nseq);
void align(hconf_t *conf, hctx_t *ctx, readinfo_t *readinfo, void *idx, runconf_t *runconf) {
    char *input = hconf_input(conf);
    const int I = hctx_mynode(ctx);
    const int N = hctx_num_nodes(ctx);

    hfile_t infile;
    hfile_init_input_stream(infile, "%s", input);
    LOG_INFO("Node %d: Aligning reads from %d input blocks.", I, ri_nblocks(*readinfo));

    hfile_t outfile;
    hfile_init_output_stream(outfile, "%s/aln-%03d.bin", hconf_tempdir(conf), I);

    const int NTHREADS = runconf->NTHREADS;
    mm_tbuf_t tbuf[NTHREADS];
    for (int i = 0; i < NTHREADS; ++i) mm_tbuf_init(&tbuf[i]);

    const size_t max_nbases_per_block = hconf_max_bases_per_batch(conf);
    size_t bufsz = max_nbases_per_block + (1 MiB);
    char *buffer = hmalloc(char, bufsz);
    shared_data_t sd = {
            .w = hconf_window_len(conf), .k = hconf_kmer_len(conf),
            .readinfo = readinfo, .buffer = buffer, .idx = idx, .tbuf = tbuf
    };
    void *tpool = runconf->tpool;
    size_t n_aln = 0, n_hits = 0;

    clist_t *list = clist_create();
    u128_v blocksizes = {0, 0, 0};

    const int my_lst_blk = ri_blkoffset(*readinfo) + ri_nblocks(*readinfo);
    int max_blocks_per_node = readinfo->blockv.n / N + (readinfo->blockv.n % N > 0);
    int round = 0;
    for (int i = ri_blkoffset(*readinfo); i < my_lst_blk; ++round, ++i) {
        sd.block = ri_block(*readinfo, i);
        uint32_t starting_rid = ri_block_startid(*readinfo, i),
                nseq = ri_block_nreads(*readinfo, i),
                ending_rid = starting_rid + nseq - 1;
        long fileoffset = ri_readoffset(*readinfo, starting_rid);
        size_t readsize = ri_readoffset(*readinfo, ending_rid) +
                          ri_readlen(*readinfo, ending_rid) + 1 - fileoffset;
        LOG_INFO("Node %d: Processing block %d of %zd with %u reads (rid %u to %u) and %zd bytes..",
                 I, i+1, readinfo->blockv.n, nseq, starting_rid, ending_rid, readsize);
        if (bufsz < readsize) {
            buffer = hrealloc(buffer, char, readsize);
            bufsz = readsize;
        }
        if (!hfile_reset(infile, fileoffset))
            LOG_ERR("Node %d: Failed to seek %ld bytes.", I, fileoffset);
        if (readsize != hfile_read(infile, buffer, char, readsize))
            LOG_ERR("Node %d: Failed to read %zd bytes from offset %ld.", I, readsize, fileoffset);
        sd.buffer = buffer;
        if (N > 1) process_block(&sd, runconf, nseq);
        else kt_forpool(runconf->tpool, smp_aligner_fn, &sd, nseq);
        /// collate hits
        for (int t = 0; t < NTHREADS; ++t) {
            if (sd.tbuf[t].hits.n) {
                clist_push(list, ma_hit_t, sd.tbuf[t].hits.a, sd.tbuf[t].hits.n);
                n_hits += sd.tbuf[t].hits.n;
                sd.tbuf[t].hits.n = 0;
                while (clist_length(list) > 1) {
                    u_char *data; size_t nbytes, count;
                    clist_pop_compressed(list, &data, &nbytes, &count);
                    hfile_write(outfile, data, u_char, nbytes);
                    hfree(data);
                    u128_t b = {nbytes, count};
                    kv_push(u128_t, blocksizes, b);
                }
            }
        }
    }
    for (; round < max_blocks_per_node; ++round) {
        LOG_INFO("Node %d: Waiting..", I);
        anonymous_barrier;
        idx_wait(idx);
        anonymous_barrier;
    }
    for (int t = 0; t < NTHREADS; ++t) {
        n_aln += sd.tbuf[t].n_aln;
        mm_tbuf_clear(&tbuf[t]);
    }
    hfile_close_stream(outfile);
    hfile_init_output_stream(outfile, "%s/alninfo-%03d.bin", hconf_tempdir(conf), I);
    hfile_write(outfile, &blocksizes.n, size_t, 1);
    hfile_write(outfile, blocksizes.a, u128_t, blocksizes.n);
    hfile_close_stream(outfile);
    LOG_INFO("Node %d: Generated %zd hits from %zd alignments.", I, n_hits, n_aln);
    hfree(buffer);

    /// free index memory
    idx_destroy(&idx);
}

extern void mm_sketch(const char *str, int len, int w, int k, uint64_t goff, u128_v *p);
#include "utils.h"
#define RADIUS 500
#define MAXGAP 10000
#define MINCNT 4
#define MERGEFRAC 0
#define MINMATCH 40

#include "ksort.h"
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, u128_t, sort_key_128x, 8)
#define sort_key_64(a) (a)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)
#define lt_low32(a, b) ((uint32_t)(a) < (uint32_t)(b))
KSORT_INIT(low32lt, uint64_t, lt_low32)
#define gt_low32(a, b) ((uint32_t)(a) > (uint32_t)(b))
KSORT_INIT(low32gt, uint64_t, gt_low32)

static void sketcher_fn(void *s, long i, int tid);
static void aligner_fn(void *s, long i, int tid);
static void req_setup_fn(void *s, long ignore, int tid);
static void get_reg(mm_tbuf_t *b, u128_t *mm, int K);
static inline void push_intv(u128_v *intv, int start, int end, float merge_frac);
static void proc_intv(mm_tbuf_t *b, u128_t *mm, int which, int k, int min_cnt, int max_gap);
static void save_hits(const mm_reg1_t *regs, int n_regs, mm_tbuf_t *b, readinfo_t *readinfo, uint my_rid, uint my_rdlen);

static void smp_aligner_fn(void *s, long i, int tid) {
    shared_data_t *sd = (shared_data_t*)s;
    readinfo_t *readinfo = sd->readinfo;
    mm_tbuf_t *b = &sd->tbuf[tid];
    u128_v *mmv = &sd->tbuf[tid].v_mm;
    uint64_t block = sd->block;
    uint32_t starting_rid = (uint32_t)(block>>32u);
    uint32_t my_rid = starting_rid + i;
    uint32_t my_rdlen = ri_readlen(*readinfo, my_rid);
    size_t off = ri_readoffset(*readinfo, my_rid) - ri_readoffset(*readinfo, starting_rid);
    char *read = sd->buffer + off;
    mmv->n = 0;
    mm_sketch(read, my_rdlen, sd->w, sd->k, 0, mmv);
    if (mmv->n) {
        for (int j = 0; j < mmv->n; ++j)
            mmv->a[j].y |= ((uint64_t)my_rid)<<32u;
    }
    b->coef.n = 0;
    uint16_t n;
    const uint64_t *r;
    for (size_t j = 0; j < mmv->n; ++j) {
        uint32_t pos = mm_pos(mmv->a[j]);
        int32_t qpos = pos >> 1u,
                strand = pos & 1u;
        r = idx_get_vals(sd->idx, tid, mmv->a[j].x, &n);
        if (n) {
            for (int k = 0; k < n; ++k) {
                int32_t rpos = ((uint32_t)r[k]) >> 1u;
                uint32_t rid = (uint32_t)(r[k] >> 32u);

                /// if my read-id is the same as rid, continue (because we ignore self matches)
                if (my_rid == rid && rpos == qpos) continue;

                ///if my read-id < rid, continue (because in all-vs-all mapping we only generate unordered read pairs)
                if (hashint32(my_rid) < hashint32(rid)) continue;

                u128_t *p;
                kv_pushp(u128_t, b->coef, &p);

                if ((r[k]&1) == strand) { /// forward strand
                    p->x = (uint64_t)r[k] >> 32 << 32 | (0x80000000U + rpos - qpos);
                    p->y = (uint64_t)j << 32 | rpos;
                } else { /// reverse strand
                    p->x = (uint64_t)r[k] >> 32 << 32 | (rpos + qpos) | 1ULL<<63;
                    p->y = (uint64_t)j << 32 | rpos;
                }
            }
        }
        mmv->a[j].y = mmv->a[j].y << 32 >> 32; /// clear the rid field
    }
    radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);
    kv_clear(b->reg);
    get_reg(b, mmv->a, sd->k);
    uint32_t qlen = ri_readlen(*readinfo, my_rid);
    save_hits(b->reg.a, b->reg.n, b, readinfo, my_rid, qlen);
}
static void process_block(shared_data_t *sd, runconf_t *runconf, size_t nseq) {
    /// sketch
    for (int i = 0; i < runconf->NTHREADS; ++i)
        sd->tbuf[i].v_mm.n = 0;
    kt_forpool(runconf->tpool, sketcher_fn, sd, nseq);

    /// setup
    kt_forpool(runconf->tpool, req_setup_fn, sd, runconf->NTHREADS);

    /// request values from peers
    anonymous_barrier;
    idx_send_requests(sd->idx);
    anonymous_barrier;

    /// align
    kt_forpool(runconf->tpool, aligner_fn, sd, runconf->NTHREADS);
}

static void sketcher_fn(void *s, long i, int tid) {
    shared_data_t *sd = (shared_data_t*)s;
    readinfo_t *readinfo = sd->readinfo;
    u128_v *mmv = &sd->tbuf[tid].v_mm;
    uint64_t block = sd->block;
    uint32_t starting_rid = (uint32_t)(block>>32u);
    uint32_t my_rid = starting_rid + i;
    uint32_t my_rdlen = ri_readlen(*readinfo, my_rid);
    size_t off = ri_readoffset(*readinfo, my_rid) - ri_readoffset(*readinfo, starting_rid);
    char *read = sd->buffer + off;
    size_t start = mmv->n;
    mm_sketch(read, my_rdlen, sd->w, sd->k, 0, mmv);
    size_t count = mmv->n - start;
    if (count) {
        for (int j = start; j < mmv->n; ++j)
            mmv->a[j].y |= (((uint64_t)my_rid)<<32u);
    }
}
static void req_setup_fn(void *s, long i, int ignore) {
    shared_data_t *sd = (shared_data_t*)s;
    u128_v *mmv = &sd->tbuf[i].v_mm;
    idx_setup_request(sd->idx, i, mmv->a, mmv->n);
}
static void aligner_fn(void *s, long c_id, int tid) {
    shared_data_t *sd = (shared_data_t*)s;
    idx_setup_response(sd->idx, c_id);
    uint n_reads = idx_nreads(sd->idx, c_id);
    for (int i = 0; i < n_reads; ++i) {
        uint nmm, qid;
        u128_t *mm = idx_get_mm(sd->idx, c_id, i, &nmm, &qid);
#ifdef SANITY_CHECKS
        uint lastpos = mm_pos(mm[0]);
        assert(mm_rid(mm[0]) == qid);
        for (int j = 1; j < nmm; ++j) {
            assert(mm_rid(mm[j]) == qid);
            assert(mm_pos(mm[j]) >= lastpos);
            lastpos = mm_pos(mm[j]);
        }
#endif
        mm_tbuf_t *b = &sd->tbuf[tid];
        b->coef.n = 0;
        uint16_t n;
        const uint64_t *r;
        for (size_t j = 0; j < nmm; ++j) {
            uint32_t pos = mm_pos(mm[j]);
            int32_t qpos = pos >> 1u,
                    strand = pos & 1u;
            r = idx_get_vals(sd->idx, c_id, mm[j].x, &n);
            if (n) {
                for (int k = 0; k < n; ++k) {
                    int32_t rpos = ((uint32_t)r[k]) >> 1u;
                    uint32_t rid = (uint32_t)(r[k] >> 32u);

                    /// if my read-id is the same as rid, continue (because we ignore self matches)
                    if (qid == rid && rpos == qpos) continue;

                    ///if my read-id > rid, continue (because in all-vs-all mapping we only generate unordered read pairs)
                    if (hashint32(qid) < hashint32(rid)) continue;

                    u128_t *p;
                    kv_pushp(u128_t, b->coef, &p);

                    if ((r[k]&1) == strand) { /// forward strand
                        p->x = (uint64_t)r[k] >> 32 << 32 | (0x80000000U + rpos - qpos);
                        p->y = (uint64_t)j << 32 | rpos;
                    } else { /// reverse strand
                        p->x = (uint64_t)r[k] >> 32 << 32 | (rpos + qpos) | 1ULL<<63;
                        p->y = (uint64_t)j << 32 | rpos;
                    }
                }
            }
            mm[j].y = mm[j].y << 32 >> 32; /// clear the rid field
        }
        radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);
        kv_clear(b->reg);
        get_reg(b, mm, sd->k);
        readinfo_t *readinfo = sd->readinfo;
        uint32_t qlen = ri_readlen(*readinfo, qid);
        save_hits(b->reg.a, b->reg.n, b, readinfo, qid, qlen);
    }
}
static void save_hits(const mm_reg1_t *regs, int n_regs, mm_tbuf_t *b, readinfo_t *readinfo, uint my_rid, uint my_rdlen) {
    for (int j = 0; j < n_regs; ++j) {
        const mm_reg1_t *r = &regs[j];
        if (r->len < MINMATCH) continue;
        b->n_aln++;
        paf_rec_t p = {.qid=my_rid, .ql=my_rdlen, .qs=r->qs, .qe=r->qe,
                .rev=r->rev, .tid=r->rid, .tl=ri_readlen(*readinfo, r->rid), .ts=r->rs, .te=r->re,
                .ml = r->len, .bl=(r->re - r->rs > r->qe - r->qs ? r->re - r->rs : r->qe - r->qs)};
        if (p.qe - p.qs < min_span || p.te - p.ts < min_span || p.ml < min_match) continue;
        ma_hit_t *h;
        kv_pushp(ma_hit_t, b->hits, &h);
        h->qns = (uint64_t)p.qid <<32u | p.qs;
        h->qe = p.qe;
        h->tn = p.tid;
        h->ts = p.ts, h->te = p.te, h->rev = p.rev, h->ml = p.ml, h->bl = p.bl;
        if (bi_dir && p.qid != p.tid) {
            kv_pushp(ma_hit_t, b->hits, &h);
            h->qns = (uint64_t)p.tid <<32u | p.ts;
            h->qe = p.te;
            h->tn = p.qid;
            h->ts = p.qs, h->te = p.qe, h->rev = p.rev, h->ml = p.ml, h->bl = p.bl;
        }
        assert(h->qns>>32 < ri_nreads(*readinfo));
    }
}
static void get_reg(mm_tbuf_t *b, u128_t *mm, const int K) {
    const uint64_t v_kept = ~(1ULL<<31), v_dropped = 1ULL<<31;
    u128_v *c = &b->coef;
    int i, j, start = 0, iso_dist = RADIUS * 2;

    if (kv_size(*c) < MINCNT) return;

    /// identify (possibly overlapping) intervals within _radius_; an interval is a cluster of hits
    b->intv.n = 0;
    for (i = 1; i < kv_size(*c); ++i) {
        if (c->a[i].x - c->a[start].x > RADIUS) {
            if (i - start >= MINCNT) push_intv(&b->intv, start, i, MERGEFRAC);
            for (++start; start < i && c->a[i].x - c->a[start].x > RADIUS; ++start);
        }
    }
    if (i - start >= MINCNT) push_intv(&b->intv, start, i, MERGEFRAC);

    /// sort by the size of the interval
    radix_sort_128x(b->intv.a, b->intv.a + b->intv.n);

    /// generate hits, starting from the largest interval
    b->reg2mini.n = 0;
    for (i = b->intv.n - 1; i >= 0; --i) proc_intv(b, mm, i, K, MINCNT, MAXGAP);
}
/// merge or add a Hough interval; only used by get_reg()
static inline void push_intv(u128_v *intv, int start, int end, float merge_frac) {
    u128_t *p;
    if (intv->n > 0) { /// test overlap
        int last_start, last_end, min;
        p = &intv->a[intv->n-1];
        last_start = p->y, last_end = p->x + last_start;
        min = end - start < last_end - last_start? end - start : last_end - last_start;
        if (last_end > start && last_end - start > min * merge_frac) { /// large overlap; then merge
            p->x = end - last_start;
            return;
        }
    }
    kv_pushp(u128_t, *intv, &p); /// a new interval
    p->x = end - start, p->y = start;
}
static void proc_intv(mm_tbuf_t *b, u128_t *mm, int which, int k, int min_cnt, int max_gap) {
    int i, j, l_lis, rid = -1, rev = 0, start = b->intv.a[which].y, end = start + b->intv.a[which].x;

    /// make room for arrays needed by LIS (longest increasing sequence)
    if (end - start > b->m) {
        b->m = end - start;
        kv_roundup32(b->m);
        b->a = hrealloc(b->a, uint64_t, b->m);
        b->b = hrealloc(b->b, size_t, b->m);
        b->p = hrealloc(b->p, size_t, b->m);
    }

    /// prepare the input array _a_ for LIS
    b->n = 0;
    for (i = start; i < end; ++i)
        if (b->coef.a[i].x != UINT64_MAX)
            b->a[b->n++] = b->coef.a[i].y, rid = b->coef.a[i].x << 1 >> 33, rev = b->coef.a[i].x >> 63;
    if (b->n < min_cnt) return;
    radix_sort_64(b->a, b->a + b->n);

    /// find the longest increasing sequence
    l_lis = rev? ks_lis_low32gt(b->n, b->a, b->b, b->p) : ks_lis_low32lt(b->n, b->a, b->b, b->p); /// LIS
    if (l_lis < min_cnt) return;
    for (i = 1, j = 1; i < l_lis; ++i) /// squeeze out minimizers reused in the LIS sequence
        if (b->a[b->b[i]]>>32 != b->a[b->b[i-1]]>>32)
            b->a[b->b[j++]] = b->a[b->b[i]];
    l_lis = j;
    if (l_lis < min_cnt) return;

    /// convert LISes to regions; possibly break an LIS at a long gaps
    for (i = 1, start = 0; i <= l_lis; ++i) {
        int32_t qgap = i == l_lis? 0 :
                       ((uint32_t)mm[b->a[b->b[i]]>>32].y>>1) - ((uint32_t)mm[b->a[b->b[i-1]]>>32].y>>1);
        if (i == l_lis || (qgap > max_gap && abs((int32_t)b->a[b->b[i]] - (int32_t)b->a[b->b[i-1]]) > max_gap)) {
            if (i - start >= min_cnt) {
                uint32_t lq = 0, lr = 0, eq = 0, er = 0, sq = 0, sr = 0;
                mm_reg1_t *r;
                kv_pushp(mm_reg1_t, b->reg, &r);
                r->rid = rid, r->rev = rev, r->cnt = i - start, r->rep = 0;
                r->qs = ((uint32_t)mm[b->a[b->b[start]]>>32].y>>1) - (k - 1);
                r->qe = ((uint32_t)mm[b->a[b->b[i-1]]>>32].y>>1) + 1;
                r->rs = rev? (uint32_t)b->a[b->b[i-1]] : (uint32_t)b->a[b->b[start]];
                r->re = rev? (uint32_t)b->a[b->b[start]] : (uint32_t)b->a[b->b[i-1]];
                r->rs -= k - 1;
                r->re += 1;
                for (j = start; j < i; ++j) { /// count the number of times each minimizer is used
                    int jj = b->a[b->b[j]]>>32;
                    mm[jj].y += 1ULL<<32;
                    kv_push(uint32_t, b->reg2mini, jj); /// keep minimizer<=>reg mapping for derep
                }
                for (j = start; j < i; ++j) { /// compute ->len
                    uint32_t q = ((uint32_t)mm[b->a[b->b[j]]>>32].y>>1) - (k - 1);
                    uint32_t r = (uint32_t)b->a[b->b[j]];
                    r = !rev? r - (k - 1) : (0x80000000U - r);
                    if (r > er) lr += er - sr, sr = r, er = sr + k;
                    else er = r + k;
                    if (q > eq) lq += eq - sq, sq = q, eq = sq + k;
                    else eq = q + k;
                }
                lr += er - sr, lq += eq - sq;
                r->len = lr < lq? lr : lq;
            }
            start = i;
        }
    }
}