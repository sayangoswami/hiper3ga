//
// Created by sayan on 9/13/19.
//

#include "hiper3ga.h"
#include "kthread.h"
#include "hfile.h"
#include "clist.h"
#include "metadata.h"

typedef struct {
    int w, k;
    readinfo_t *readinfo;
    uint64_t block;
    char *buffer;
    u128_v *mmv;
    clist_t **lists;
} shared_data_t;

static void sketcher_fn(void *s, long i, int tid);
void mm_sketch(const char *str, int len, int w, int k, uint64_t goff, u128_v *p);

/**
 * Takes a set of reads and emits <f,r,p> tuples where f = minimizer, r = read-id, p = position of f in r
 * @param conf
 * @param split
 * @param ctx
 */
clist_t* mapper(hconf_t *conf, hctx_t *ctx, readinfo_t *readinfo, runconf_t *runconf) {
    char *input = hconf_input(conf);
    const int I = hctx_mynode(ctx);

    hfile_t infile;
    hfile_init_input_stream(infile, "%s", input);
    LOG_INFO("Node %d: Generating fingerprints from %d input blocks.", I, ri_nblocks(*readinfo));
    uint my_total_nseq = 0;
    const size_t max_nbases_per_block = hconf_max_bases_per_batch(conf);
    size_t bufsz = max_nbases_per_block + (1 MiB);
    char *buffer = hmalloc(char, bufsz);
    const int NTHREADS = runconf->NTHREADS;
    u128_v *mmv = hmalloc(u128_v, NTHREADS);
    clist_t *lists[NTHREADS];
    for (int i = 0; i < NTHREADS; ++i) lists[i] = clist_create();
    shared_data_t sd = {
                .readinfo = readinfo, .mmv = mmv, .lists = lists,
                .w = hconf_window_len(conf), .k = hconf_kmer_len(conf),
    };
    void *tpool = runconf->tpool;

    const int my_lst_blk = ri_blkoffset(*readinfo) + ri_nblocks(*readinfo);
    for (int i = ri_blkoffset(*readinfo); i < my_lst_blk; ++i) {
        sd.block = ri_block(*readinfo, i);
        uint32_t starting_rid = ri_block_startid(*readinfo, i),
                 nseq = ri_block_nreads(*readinfo, i),
                 ending_rid = starting_rid + nseq-1;
        long fileoffset = ri_readoffset(*readinfo, starting_rid);
        size_t readsize = ri_readoffset(*readinfo, ending_rid) +
                ri_readlen(*readinfo, ending_rid) + 1 - fileoffset;
        LOG_INFO("Node %d: Processing block %d with %u reads (rid %u to %u) and %zd bytes..",
                I, i, nseq, starting_rid, ending_rid, readsize);
        if (bufsz < readsize) {
            buffer = hrealloc(buffer, char, readsize);
            bufsz = readsize;
        }
        if (!hfile_reset(infile, fileoffset))
            LOG_ERR("Node %d: Failed to seek %ld bytes.", I, fileoffset);
        if (readsize != hfile_read(infile, buffer, char, readsize))
            LOG_ERR("Node %d: Failed to read %zd bytes from offset %ld.", I, readsize, fileoffset);
        sd.buffer = buffer;
        kt_forpool(tpool, sketcher_fn, &sd, nseq);
        my_total_nseq += nseq;
    }
    hfile_close_stream(infile); hfree(buffer);
    clist_t *out = clist_create();
    for (int i = 0; i < NTHREADS; ++i) {
        clist_append(out, lists[i]);
        kv_destroy(sd.mmv[i]);
    }
    hfree(mmv);

    LOG_INFO("Node %d: From %u reads, generated %zd tuples (%.2fx compressed) in %zd blocks.",
             I, my_total_nseq, clist_count(out), clist_compfactor(out, u128_t), clist_length(out));
    return out;
}

static void sketcher_fn(void *s, long i, int tid) {
    shared_data_t *sd = (shared_data_t*)s;
    readinfo_t *readinfo = sd->readinfo;
    u128_v *mmv = &sd->mmv[tid];
    uint64_t block = sd->block;
    uint32_t starting_rid = (uint32_t)(block>>32u);
    uint32_t my_rid = starting_rid + i;
    size_t off = ri_readoffset(*readinfo, my_rid) - ri_readoffset(*readinfo, starting_rid);
    char *read = sd->buffer + off;
    mmv->n = 0;
    mm_sketch(read, ri_readlen(*readinfo, my_rid), sd->w, sd->k, readinfo->cumulen[my_rid], mmv);
    if (mmv->n)
        clist_push(sd->lists[tid], u128_t, mmv->a, mmv->n);
}

unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param goff   global offset of the read in the string obtained by appending all reads
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = (goff + lastPos)<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(const char *str, int len, int w, int k, uint64_t goff, u128_v *p) {
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
    int i, j, l, buf_pos, min_pos;
    u128_t *buf, min = {UINT64_MAX, UINT64_MAX };

    assert(len > 0 && w > 0 && k > 0);
    buf = (u128_t*)alloca(w * 16);
    memset(buf, 0xff, w * 16);

    for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)str[i]];
        u128_t info = {UINT64_MAX, UINT64_MAX };
        if (c < 4) { /// not an ambiguous base
            int z;
            kmer[0] = (kmer[0] << 2 | c) & mask;           /// forward k-mer
            kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; /// reverse k-mer
            if (kmer[0] == kmer[1]) continue; /// skip "symmetric k-mers" as we don't know it strand
            z = kmer[0] < kmer[1]? 0 : 1; /// strand
            if (++l >= k)
                info.x = hash64(kmer[z], mask), info.y = (goff+i)<<1 | z;
        } else l = 0;
        buf[buf_pos] = info; /// need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
        if (l == w + k - 1) { /// special case for the first window - because identical k-mers are not stored yet
            for (j = buf_pos + 1; j < w; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) kv_push(u128_t, *p, buf[j]);
            for (j = 0; j < buf_pos; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) kv_push(u128_t, *p, buf[j]);
        }
        if (info.x <= min.x) { /// a new minimum; then write the old min
            if (l >= w + k) kv_push(u128_t, *p, min);
            min = info, min_pos = buf_pos;
        } else if (buf_pos == min_pos) { /// old min has moved outside the window
            if (l >= w + k - 1) kv_push(u128_t, *p, min);
            for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) /// the two loops are necessary when there are identical k-mers
                if (min.x >= buf[j].x) min = buf[j], min_pos = j; /// >= is important s.t. min is always the closest k-mer
            for (j = 0; j <= buf_pos; ++j)
                if (min.x >= buf[j].x) min = buf[j], min_pos = j;
            if (l >= w + k - 1) { /// write identical k-mers
                for (j = buf_pos + 1; j < w; ++j) /// these two loops make sure the output is sorted
                    if (min.x == buf[j].x && min.y != buf[j].y) kv_push(u128_t, *p, buf[j]);
                for (j = 0; j <= buf_pos; ++j)
                    if (min.x == buf[j].x && min.y != buf[j].y) kv_push(u128_t, *p, buf[j]);
            }
        }
        if (++buf_pos == w) buf_pos = 0;
    }
    if (min.x != UINT64_MAX)
        kv_push(u128_t, *p, min);
}