//
// Created by sayan on 6/8/20.
//

#include "hiper3ga.h"
#include "ksort.h"
#include "metadata.h"
#include "kthread.h"
#include "clist.h"

typedef struct {
    int qid;
    ma_sub_t sub;
} qsub_t;
typedef kvec_t(qsub_t) qsub_v;

#define ma_hit_key(a) ((a).qns)
KRADIX_SORT_INIT(hit, ma_hit_t, ma_hit_key, 8)
KSORT_INIT_GENERIC(uint32_t)

#define CRUDE_SUB 1
#define CRUDE_CUT 2
#define FINE_CUT 3
#define MARK_CONTAINED 4
#define MARK_UNUSED 5
#define REMOVE_CONTAINED 6

size_t ma_hit_sub(uint32_t n_sub, qsub_v *qsub, uint64_t n, const ma_hit_t *a, int end_clip);
size_t ma_hit_cut(const ma_sub_t *sub, uint64_t n, ma_hit_t *a);
size_t ma_hit_flt(const ma_sub_t *sub, uint64_t n, ma_hit_t *a, double *cov);
static void collate_sub(hctx_t *ctx, qsub_v *qsub);
void ma_sub_merge(size_t n_sub, ma_sub_t *a, const ma_sub_t *b);
void ma_hit_mark_deleted(const ma_sub_t *sub, size_t n, const ma_hit_t *a, hbitvec_t *bv);
void ma_hit_mark_unused(size_t n, const ma_hit_t *a, hbitvec_t *bv);
size_t ma_hit_remove_contained(const int *map, size_t n, ma_hit_t *a);
static void collate_deleted(hctx_t *ctx, hbitvec_t *bv, ma_sub_t *sub, sd_seq_t *seq);
static int* collate_unused(hctx_t *ctx, hbitvec_t *bv, sd_seq_t *seq, uint32_t *new_n_seq_p);
static void readselector_fn(void *s, long id, int tid);

typedef struct {
    int step;
    uint32_t nseq;
    qsub_v qsub[ASM_NPTN];
    clist_t **lists;
    size_t n_remained[ASM_NPTN];
    ma_sub_t *sub;
    ma_sub_t *sub2;
    size_t n_hits[ASM_NPTN];
    hbitvec_t *bv;
    int *map;
    ma_hit_v hitbuf[ASM_NPTN];
} shared_data_t;

uint32_t
readselector(clist_t **lists, sd_seq_t *seq, uint32_t nseq, hctx_t *ctx, runconf_t *runconf, ma_sub_t **sub_p) {
    shared_data_t sd = {.lists = lists, .nseq = nseq};
    for (int i = 0; i < ASM_NPTN; ++i) {
        kv_init(sd.qsub[i]); sd.n_remained[i] = 0; sd.n_hits[i] = 0; kv_init(sd.hitbuf[i]);
    }

    /// crude sub
    sd.step = CRUDE_SUB;
    kt_forpool(runconf->tpool, readselector_fn, &sd, ASM_NPTN);
    collate_sub(ctx, sd.qsub);

    /// crude cut
    sd.step = CRUDE_CUT;
    sd.sub = hmalloc(ma_sub_t, nseq);
    for (int i = 0; i < sd.qsub[0].n; ++i) sd.sub[sd.qsub[0].a[i].qid] = sd.qsub[0].a[i].sub;
    kv_destroy(sd.qsub[0]);
    kt_forpool(runconf->tpool, readselector_fn, &sd, ASM_NPTN);

    /// fine sub
    collate_sub(ctx, sd.qsub);

    /// fine cut
    sd.step = FINE_CUT;
    sd.sub2 = hmalloc(ma_sub_t, nseq);
    for (int i = 0; i < sd.qsub[0].n; ++i) sd.sub2[sd.qsub[0].a[i].qid] = sd.qsub[0].a[i].sub;
    kv_destroy(sd.qsub[0]);
    kt_forpool(runconf->tpool, readselector_fn, &sd, ASM_NPTN);

    /// ma_sub_merge
    ma_sub_merge(nseq, sd.sub, sd.sub2); // TODO can we (or do we need to) distribute this?
    hfree(sd.sub2);

    /// ma_hit_contained - mark deleted/unused
    sd.bv = hbitvec_new(nseq);
    sd.step = MARK_CONTAINED;
    kt_forpool(runconf->tpool, readselector_fn, &sd, ASM_NPTN);
    collate_deleted(ctx, sd.bv, sd.sub, seq);

    hbitvec_clear_all(sd.bv);
    sd.step = MARK_UNUSED;
    kt_forpool(runconf->tpool, readselector_fn, &sd, ASM_NPTN);
    uint32_t new_nseq;
    sd.map = collate_unused(ctx, sd.bv, seq, &new_nseq);

    sd.step = REMOVE_CONTAINED;
    for (int i = 0; i < nseq; ++i)
        if (sd.map[i] >= 0) sd.sub[sd.map[i]] = sd.sub[i];
    kt_forpool(runconf->tpool, readselector_fn, &sd, ASM_NPTN);

    size_t total_hits = 0;
    for (int i = 0; i < ASM_NPTN; ++i)
        total_hits += clist_count(lists[i]);
    hfree(sd.map); hbitvec_destroy(&sd.bv);
    for (int i = 0; i < ASM_NPTN; ++i) kv_destroy(sd.hitbuf[i]);
    LOG_INFO("%d sequences and %zd hits remain after containment removal.", new_nseq, total_hits);

    ////////////////////////////////// @@@ DONE @@@ //////////////////////////////////
    *sub_p = sd.sub;
    return new_nseq;
}

static void readselector_fn(void *s, long id, int ignore) {
    shared_data_t *sd = (shared_data_t*)s;
    uint32_t n_seq = sd->nseq;
    qsub_v *qsub = &sd->qsub[id];
    ma_hit_v *hitbuf = &sd->hitbuf[id];

    if (sd->step == CRUDE_SUB) {
        hitbuf->n = clist_count(sd->lists[id]); /** number of hits */
        kv_reserve(ma_hit_t, *hitbuf, hitbuf->n);   /// array of hits
        long offset = 0, count;
        while ((count = clist_pop(sd->lists[id], ma_hit_t, hitbuf->a + offset, hitbuf->n - offset)) != 0)
            offset += count;
        clist_destroy(sd->lists[id]);
        assert(offset == hitbuf->n);

        /// sort hits
        radix_sort_hit(hitbuf->a, hitbuf->a + hitbuf->n);
        ////////////////// crude ma_hit_sub() //////////////////
        sd->n_remained[id] = ma_hit_sub(n_seq, qsub, hitbuf->n, hitbuf->a, 0);

    } else if (sd->step == CRUDE_CUT) {

        ma_sub_t *sub = sd->sub;

        ////////////////// crude ma_hit_cut() //////////////////
        size_t n_hits = ma_hit_cut(sub, hitbuf->n, hitbuf->a);

        ////////////////// crude ma_hit_flt() //////////////////
        double cov;
        n_hits = ma_hit_flt(sub, n_hits, hitbuf->a, &cov);

        ////////////////// fine ma_hit_sub() //////////////////
        sd->n_remained[id] = ma_hit_sub(n_seq, qsub, n_hits, hitbuf->a, min_span/2);

        sd->n_hits[id] = n_hits;

    } else if (sd->step == FINE_CUT) {

        ////////////////// fine ma_hit_cut() //////////////////
        ma_sub_t *sub2 = sd->sub2;
        sd->n_hits[id] = ma_hit_cut(sub2, sd->n_hits[id], hitbuf->a);

    } else if (sd->step == MARK_CONTAINED) {

        ////////////////// mark deleted //////////////////
        ma_hit_mark_deleted(sd->sub, sd->n_hits[id], hitbuf->a, sd->bv);

    } else if (sd->step == MARK_UNUSED) {

        ////////////////// mark unused //////////////////
        ma_hit_mark_unused(sd->n_hits[id], hitbuf->a, sd->bv);

    } else if (sd->step == REMOVE_CONTAINED) {

        size_t n_remained = ma_hit_remove_contained(sd->map, sd->n_hits[id], hitbuf->a);
        kv_reserve(ma_hit_t, *hitbuf, 2*n_remained);
        sd->lists[id] = clist_create();
        clist_push(sd->lists[id], ma_hit_t, hitbuf->a, n_remained);
    }
}

/**
 * gather from workers, and send to master if N > 1
 * @param globe
 * @param qsub
 */
static void collate_sub(hctx_t *ctx, qsub_v *qsub) {
    const int N = hctx_num_nodes(ctx);
    const int I = hctx_mynode(ctx);
    size_t count = 0;
    for (int i = 0; i < ASM_NPTN; ++i) count += kv_size(qsub[i]);
    kv_reserve(qsub_t, qsub[0], count);
    for (size_t i = 1, off = qsub[0].n; i < ASM_NPTN; ++i){
        memcpy(qsub[0].a + off, qsub[i].a, qsub[i].n * sizeof(qsub_t));
        off += qsub[i].n;
        qsub[0].n = off;
        kv_destroy(qsub[i]);
    }
    if (N > 1){
        anonymous_barrier;

        /// send to master and receive from master
        size_t *nqsub_p = hctx_segaddr(ctx, I);
        qsub_t *qsubarr = hctx_segaddr(ctx, I) + sizeof(size_t);
        *nqsub_p = qsub[0].n;
        memcpy(qsubarr, qsub[0].a, qsub[0].n * sizeof(qsub_t));
        anonymous_barrier;

        if (I == 0){
            /// I am the master
            size_t rcounts[N], total_count = 0;
            gasnet_valget_handle_t cnthandles[N];

            for (uint i = 0; i < N; ++i) cnthandles[i] = gasnet_get_nb_val(i, hctx_segaddr(ctx, i), sizeof(size_t));
            for (uint i = 0; i < N; ++i) {
                rcounts[i] = gasnet_wait_syncnb_valget(cnthandles[i]);
                total_count += rcounts[i];
            }
            kv_resize(qsub_t, qsub[0], total_count);
            size_t off = 0;
            for (uint i = 0; i < N; ++i) {
                gasnet_get_nbi(qsub[0].a + off, i, hctx_segaddr(ctx, i) + sizeof(size_t), rcounts[i] * sizeof(qsub_t));
                off += rcounts[i];
            } gasnet_wait_syncnbi_gets();
            assert(off == total_count);
            qsub[0].n = total_count;
            for (uint i = 0; i < N; ++i) {
                gasnet_put_nbi_val(i, hctx_segaddr(ctx, i), total_count, sizeof(size_t));
                gasnet_put_nbi(i, hctx_segaddr(ctx, i) + sizeof(size_t), qsub[0].a, qsub[0].n * sizeof(qsub_t));
            } gasnet_wait_syncnbi_puts();

            anonymous_barrier;
        } else {
            anonymous_barrier;
            kv_resize(qsub_t, qsub[0], *nqsub_p);
            qsub[0].n = *nqsub_p;
            memcpy(qsub[0].a, qsubarr, qsub[0].n * sizeof(qsub_t));
        }
    }
}

static void collate_deleted(hctx_t *ctx, hbitvec_t *bv, ma_sub_t *sub, sd_seq_t *seq) {
    const int N = hctx_num_nodes(ctx);
    const int I = hctx_mynode(ctx);
    size_t count = 0;
    uint32_t *array = hbitvec_raw(bv, &count);
    const size_t array_size = count/32 + (count%32 != 0);
    if (N > 1) {
        size_t *count_p = hctx_segaddr(ctx, I);
        uint32_t *myarray = hctx_segaddr(ctx, I) + sizeof(size_t);
        *count_p = array_size;
        memcpy(myarray, array, array_size * sizeof(uint32_t));
        anonymous_barrier;

        if (I == 0) { /// I am the master
            uint32_t *data[N];
            size_t rcounts[N];
            gasnet_valget_handle_t cnthandles[N];
            for (uint i = 0; i < N; ++i) cnthandles[i] = gasnet_get_nb_val(i, hctx_segaddr(ctx, i), sizeof(size_t));
            for (uint i = 0; i < N; ++i) {
                rcounts[i] = gasnet_wait_syncnb_valget(cnthandles[i]);
                assert(rcounts[i] == array_size);
                data[i] = hmalloc(uint32_t, rcounts[i]);
            }
            for (uint i = 0; i < N; ++i)
                gasnet_get_nbi(data[i], i, hctx_segaddr(ctx, i) + sizeof(size_t), rcounts[i] * sizeof(uint32_t));
            gasnet_wait_syncnbi_gets();

            for (int i = 0; i < N; ++i) {
                for (size_t j = 0; j < array_size; ++j) array[j] |= data[i][j];
                hfree(data[i]);
            }

            for (uint i = 0; i < N; ++i) {
                gasnet_put_nbi_val(i, hctx_segaddr(ctx, i), array_size, sizeof(size_t));
                gasnet_put_nbi(i, hctx_segaddr(ctx, i) + sizeof(size_t), array, array_size * sizeof(uint32_t));
            } gasnet_wait_syncnbi_puts();

            anonymous_barrier;
        } else {
            anonymous_barrier;
            assert(*count_p == array_size);
            memcpy(array, myarray, array_size * sizeof(uint32_t));
        }
    }
    for (size_t i = 0; i < count; ++i) {
        size_t offset = i / 32;
        size_t bitpos = i % 32;
        if (array[offset] & (0x80000000 >> (bitpos))) {
            sub[i].del = seq[i].del = 1;
        }
    }
}

static int* squeeze(sd_seq_t *seq, uint32_t n_seq, uint32_t *new_n_seq_p) {
    int *map = hmalloc(int, n_seq);
    uint32_t i, j;
    for (i = j = 0; i < n_seq; ++i) {
        if (seq[i].del) {
            map[i] = -1;
        } else seq[j] = seq[i], map[i] = j++;
    }
    *new_n_seq_p = j;
    return map;
}

static int* collate_unused(hctx_t *ctx, hbitvec_t *bv, sd_seq_t *seq, uint32_t *new_n_seq_p) {
    const int N = hctx_num_nodes(ctx);
    const int I = hctx_mynode(ctx);
    size_t count = 0;
    uint32_t *array = hbitvec_raw(bv, &count);
    const size_t array_size = count/32 + (count%32 != 0);
    if (N > 1) {
        /// scatter & gather
        size_t *count_p = hctx_segaddr(ctx, I);
        uint32_t *myarray = hctx_segaddr(ctx, I) + sizeof(size_t);
        *count_p = array_size;
        memcpy(myarray, array, array_size * sizeof(uint32_t));
        anonymous_barrier;
        if (I == 0) { /// I am the master
            uint32_t *data[N];
            size_t rcounts[N];
            gasnet_valget_handle_t cnthandles[N];
            for (uint i = 0; i < N; ++i) cnthandles[i] = gasnet_get_nb_val(i, hctx_segaddr(ctx, i), sizeof(size_t));
            for (uint i = 0; i < N; ++i) {
                rcounts[i] = gasnet_wait_syncnb_valget(cnthandles[i]);
                assert(rcounts[i] == array_size);
                data[i] = hmalloc(uint32_t, rcounts[i]);
            }
            for (uint i = 0; i < N; ++i)
                gasnet_get_nbi(data[i], i, hctx_segaddr(ctx, i) + sizeof(size_t), rcounts[i] * sizeof(uint32_t));
            gasnet_wait_syncnbi_gets();
            for (int i = 0; i < N; ++i) {
                for (size_t j = 0; j < array_size; ++j) array[j] |= data[i][j];
                hfree(data[i]);
            }
            for (uint i = 0; i < N; ++i) {
                gasnet_put_nbi_val(i, hctx_segaddr(ctx, i), array_size, sizeof(size_t));
                gasnet_put_nbi(i, hctx_segaddr(ctx, i) + sizeof(size_t), array, array_size * sizeof(uint32_t));
            } gasnet_wait_syncnbi_puts();

            anonymous_barrier;
        } else {
            anonymous_barrier;
            assert(*count_p == array_size);
            memcpy(array, myarray, array_size * sizeof(uint32_t));
        }
    }
    for (size_t i = 0; i < count; ++i) {
        size_t offset = i / 32;
        size_t bitpos = i % 32;
        sd_seq_t *s = &seq[i];
        if (!(array[offset] & (0x80000000 >> (bitpos)))) {
            s->aux = 0;
            s->del = 1;
        }
        else s->aux = 1;
    }
    int *map = squeeze(seq, count, new_n_seq_p);
    return map;
}

void ma_sub_merge(size_t n_sub, ma_sub_t *a, const ma_sub_t *b) {
    size_t i;
    for (i = 0; i < n_sub; ++i)
        a[i].e = a[i].s + b[i].e, a[i].s += b[i].s;
}

size_t ma_hit_sub(uint32_t n_sub, qsub_v *qsub, uint64_t n, const ma_hit_t *a, int end_clip) {
    size_t i, j, last, n_remained = 0;
    kvec_t(uint32_t) b = {0,0,0};
    kv_init(*qsub);
    for (i = 1, last = 0; i <= n; ++i) {
        if (i == n || a[i].qns>>32 != a[i-1].qns>>32) { /// we come to a new query sequence
            size_t start = 0;
            int dp, qid = a[i-1].qns>>32;
            if(qid >= n_sub) LOG_ERR("qid %d is larger that nseq %u.", qid, n_sub);
            ma_sub_t max, max2;
            kv_reserve(uint32_t, b, i - last);
            b.n = 0;
            for (j = last; j < i; ++j) { /// collect all starts and ends
                uint32_t qs, qe;
                if (a[j].tn == qid || a[j].ml < a[j].bl * min_iden) continue; /// skip self match
                qs = (uint32_t)a[j].qns + end_clip, qe = a[j].qe - end_clip;
                if (qe > qs) {
                    kv_push(uint32_t, b, qs<<1);
                    kv_push(uint32_t, b, qe<<1|1);
                }
            }
            ks_introsort_uint32_t(b.n, b.a);
            max.s = max.e = max.del = max2.s = max2.e = max2.del = 0;
            for (j = 0, dp = 0; j < b.n; ++j) {
                int old_dp = dp;
                if (b.a[j]&1) --dp;
                else ++dp;
                if (old_dp < min_dp && dp >= min_dp) {
                    start = b.a[j]>>1;
                } else if (old_dp >= min_dp && dp < min_dp) {
                    int len = (b.a[j]>>1) - start;
                    if (len > max.e - max.s) max2 = max, max.s = start, max.e = b.a[j]>>1;
                    else if (len > max2.e - max2.s) max2.s = start, max2.e = b.a[j]>>1;
                }
            }
            qsub_t *qsp;
            kv_pushp(qsub_t, *qsub, &qsp);
            qsp->qid = qid;
            if (max.e - max.s > 0) {
                /// although sub is modified by multiple threads, no 2 threads modify the same location
                /*sub[qid].s = max.s - end_clip;
                sub[qid].e = max.e + end_clip;
                sub[qid].del = 0;*/
                qsp->sub.s = max.s - end_clip;
                qsp->sub.e = max.e + end_clip;
                qsp->sub.del = 0;
                ++n_remained;
            } else qsp->sub.del = 1; //sub[qid].del = 1;
            last = i;
        }
    }
    hfree(b.a);
    //assert(n_remained == kv_size(*qsub));
    return n_remained;
}

size_t ma_hit_cut(const ma_sub_t *sub, uint64_t n, ma_hit_t *a) {
    size_t i, m;
    for (i = m = 0; i < n; ++i) {
        ma_hit_t *p = &a[i];
        const ma_sub_t *rq = &sub[p->qns>>32], *rt = &sub[p->tn];
        int qs, qe, ts, te;
        if (rq->del || rt->del) continue;
        if (p->rev) {
            qs = p->te < rt->e? (uint32_t)p->qns : (uint32_t)p->qns + (p->te - rt->e);
            qe = p->ts > rt->s? p->qe : p->qe - (rt->s - p->ts);
            ts = p->qe < rq->e? p->ts : p->ts + (p->qe - rq->e);
            te = (uint32_t)p->qns > rq->s? p->te : p->te - (rq->s - (uint32_t)p->qns);
        } else {
            qs = p->ts > rt->s? (uint32_t)p->qns : (uint32_t)p->qns + (rt->s - p->ts);
            qe = p->te < rt->e? p->qe : p->qe - (p->te - rt->e);
            ts = (uint32_t)p->qns > rq->s? p->ts : p->ts + (rq->s - (uint32_t)p->qns);
            te = p->qe < rq->e? p->te : p->te - (p->qe - rq->e);
        }
        qs = (qs > rq->s? qs : rq->s) - rq->s;
        qe = (qe < rq->e? qe : rq->e) - rq->s;
        ts = (ts > rt->s? ts : rt->s) - rt->s;
        te = (te < rt->e? te : rt->e) - rt->s;
        if (qe - qs >= min_span && te - ts >= min_span) {
            p->qns = p->qns>>32<<32 | qs, p->qe = qe, p->ts = ts, p->te = te;
            a[m++] = *p;
        }
    }
    return m;
}


size_t ma_hit_flt(const ma_sub_t *sub, uint64_t n, ma_hit_t *a, double *cov) {
    size_t i, m;
    asg_arc_t t;
    uint64_t tot_dp = 0, tot_len = 0;
    for (i = m = 0; i < n; ++i) {
        ma_hit_t *h = &a[i];
        const ma_sub_t *sq = &sub[h->qns>>32], *st = &sub[h->tn];
        int r;
        if (sq->del || st->del) continue;
        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang*1.5, .5, min_ovlp*.5, &t);
        if (r >= 0 || r == MA_HT_QCONT || r == MA_HT_TCONT)
            a[m++] = *h, tot_dp += r >= 0? r : r == MA_HT_QCONT? sq->e - sq->s : st->e - st->s;
    }
    for (i = 1; i <= m; ++i)
        if (i == m || a[i].qns>>32 != a[i-1].qns>>32)
            tot_len += sub[a[i-1].qns>>32].e - sub[a[i-1].qns>>32].s;
    *cov = (double)tot_dp / tot_len;
    return m;
}

void ma_hit_mark_deleted(const ma_sub_t *sub, const size_t n, const ma_hit_t *a, hbitvec_t *bv) {
    asg_arc_t t;
    for (size_t i = 0; i < n; ++i) {
        const ma_hit_t *h = &a[i];
        const ma_sub_t *sq = &sub[h->qns>>32u], *st = &sub[h->tn];
        int r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, int_frac, min_ovlp, &t);
        if (r == MA_HT_QCONT) hbitvec_set(bv, (h->qns>>32u));
        else if (r == MA_HT_TCONT) hbitvec_set(bv, (h->tn));
    }
}

void ma_hit_mark_unused(const size_t n, const ma_hit_t *a, hbitvec_t *bv) {
    for (size_t i = 0; i < n; ++i) {
        hbitvec_set(bv, (a[i].qns>>32u));
        hbitvec_set(bv, (a[i].tn));
    }
}

size_t ma_hit_remove_contained(const int *map, size_t n, ma_hit_t *a) {
    size_t i, m;
    for (i = m = 0; i < n; ++i) {
        ma_hit_t *h = &a[i];
        int32_t qn = map[h->qns>>32], tn = map[h->tn];
        if (qn >= 0 && tn >= 0) {
            a[i].qns = (uint64_t)qn<<32 | (uint32_t)a[i].qns;
            a[i].tn = tn;
            a[m++] = a[i];
        }
    }
    return m;
}