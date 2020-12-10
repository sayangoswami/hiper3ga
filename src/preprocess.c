//
// Created by sayan on 5/14/20.
//

#include "hclasses.h"
#include "hfile.h"
#include "kthread.h"
#include "ksort.h"
#include "metadata.h"

#undef SANITY_CHECKS

#ifdef SANITY_CHECKS
char symb[] = {'*','+','#','@'};
#endif

KSORT_INIT_GENERIC(uint64_t)

typedef struct {
    u64_v *readvec;
    char *data;
    size_t file_offset, length;
    int n;
} shared_data_t;

static void reader_fn(void *s, long i, int tid);

readinfo_t * preprocess(hconf_t *conf, hctx_t *ctx) {
    char *input = hconf_input(conf);
    const int I = hctx_mynode(ctx);
    const int N = hctx_num_nodes(ctx);
    const int T = 4;
    hfile_t fin;
    hfile_init_input_stream(fin, "%s", input);
    size_t fsz = hfile_bytes(fin);
    const size_t stepsize = fsz/N;
    const size_t blksz = (I==N-1) ? (stepsize + fsz%N) : (stepsize);
    size_t start = I * stepsize, end = start + blksz;
    const size_t mapsz = 32 MiB;
    char *data = hmalloc(char, mapsz);

    u64_v *readvec = hmalloc(u64_v, T);
    shared_data_t sd = {.readvec = readvec, .n = T};
    void *tpool = kt_forpool_init(T);

    hfile_forwrd(fin, start);
    while (start < end) {
        size_t length = MIN(mapsz, end-start);
        length = hfile_read(fin, data, char, length);
        sd.data = data, sd.length = length, sd.file_offset = start;
        kt_forpool(tpool, reader_fn, &sd, T);
        start += length;
    }
    hfile_close_stream(fin);
    kt_forpool_destroy(tpool);

    /// aggregate and sort newline positions
    uint64_t *nlines_p = (uint64_t*)hctx_segaddr(ctx, I);
    uint64_t *rpos = nlines_p+1, *p = rpos;
    for (int i = 0; i < T; ++i) {
        memcpy(p, readvec[i].a, readvec[i].n * sizeof(uint64_t));
        p += readvec[i].n;
        kv_destroy(readvec[i]);
    }
    hfree(readvec);
    hfree(data);
    *nlines_p = p - rpos;
    anonymous_barrier;
    if (!I) {
        uint64_t nlinearr[N];
        for (int i = 0; i < N; ++i)
            nlinearr[i] = gasnet_get_val(i, hctx_segaddr(ctx, i), 8);
        uint64_t *dest = p;
        for (int i = 1; i < N; ++i) {
            gasnet_get_nbi(dest, i, hctx_segaddr(ctx, i) + 8, nlinearr[i] * 8);
            dest += nlinearr[i];
        }
        gasnet_wait_syncnbi_gets();
        size_t nlines = dest - rpos;
        if (nlines % 2) LOG_ERR("Expected an even number of lines.");
        ks_introsort_uint64_t(nlines, rpos);
        *nlines_p = nlines;
        for (int i = 1; i < N; ++i)
            gasnet_put_nbi(i, hctx_segaddr(ctx, i), nlines_p, (nlines+1)*8);
        gasnet_wait_syncnbi_puts();
    }
    anonymous_barrier;

    /// calculate read lengths and max read length
    size_t nreads = (*nlines_p)/2;
    u128_v *reads = hmalloc(u128_v, 1);
    kv_resize(u128_t, *reads, nreads);
    uint max_rd_len = 0;
    for (int i = 0; i < nreads; ++i) {
        reads->a[i].y = rpos[2*i+1] - rpos[2*i] - 1; // previously lengths[i]
        if (reads->a[i].y > max_rd_len)
            max_rd_len = reads->a[i].y;
    }

    /// keep every other offset (fasta file has reads every other line)
    for (int i = 0; i < nreads; ++i)
        reads->a[i].x = rpos[2*i] + 1; // previously rpos[i]

    return readinfo_create(conf, ctx, reads);
}

static void reader_fn(void *s, long i, int tid) {
    shared_data_t *sd = (shared_data_t*)s;
    u64_v *readvec = &sd->readvec[tid];
    const size_t blksz = (i==sd->n-1) ? (sd->length/sd->n + sd->length%sd->n) : (sd->length/sd->n);
    char *begin = sd->data + i*sd->length/sd->n, *end = begin + blksz;
    while (begin < end) {
        char *ret = memchr(begin, '\n', end-begin);
        if (ret) {
            kv_push(uint64_t, *readvec, (sd->file_offset + ret - sd->data));
            begin = ret+1;
        } else break;
    }
}