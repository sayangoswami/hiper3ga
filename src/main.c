//
// Created by sayan on 9/13/19.
//

#include "hiper3ga.h"
#include "kthread.h"
#include "hidx.h"
#include "metadata.h"
#include "clist.h"

extern readinfo_t *preprocess(hconf_t *conf, hctx_t *ctx);
extern clist_t* mapper(hconf_t *conf, hctx_t *ctx, readinfo_t *readinfo, runconf_t *runconf);
extern clist_t ** u128L1_partition(clist_t *listin, runconf_t *runconf);
extern clist_t ** u128L2_partition(clist_t *listin, runconf_t *runconf);
extern clist_t ** hitL1_partition(clist_t *listin, runconf_t *runconf);
extern clist_t ** hitL2_partition(clist_t *listin, runconf_t *runconf);
extern clist_t * shuffle(hctx_t *ctx, clist_t **lists);
extern void align(hconf_t *conf, hctx_t *ctx, readinfo_t *readinfo, void *idx, runconf_t *runconf);
extern uint32_t readselector(clist_t* (hits)[ASM_NPTN], sd_seq_t *seq, uint32_t nseq, hctx_t *ctx, runconf_t *runconf, ma_sub_t **sub_p);
extern void build_graph(hctx_t *ctx, const sd_seq_t *seq, uint32_t n_seq, const ma_sub_t *sub, clist_t* (hits)[ASM_NPTN], runconf_t *runconf, char *output);

int main(int argc, char **argv) {
    hctx_t *ctx = hctx_create(argc, argv);
    hconf_t *conf = hconf_create(argc, argv, ctx);

    /// create thread pool
    const int NTHREADS = hconf_nthreads(conf);
    void *tpool = kt_forpool_init(NTHREADS);
    pthread_barrier_t barrier;
    int rc = pthread_barrier_init(&barrier, NULL, NTHREADS);
    if (rc) LOG_ERR("Unable to create barrier because %s.", strerror(errno));
    runconf_t runconf = { .conf = conf, .NNODES = hctx_num_nodes(ctx),
                          .NTHREADS = NTHREADS, .barrier = &barrier, .tpool = tpool };
    snprintf(runconf.log_prefix, 64, "Node %d", hctx_mynode(ctx));

    clist_t *clist = NULL;
    clist_t **lists = NULL;

    if (streq(argv[1], "align")) {
        /// preprocess and write read-lengths to file
        readinfo_t *readinfo = preprocess(conf, ctx);
        hfile_t outfile;
        hfile_init_output_stream(outfile, "%s/rlen-%03d.bin", hconf_tempdir(conf), hctx_mynode(ctx));
        hfile_write(outfile, &ri_nreads(*readinfo), size_t, 1);
        for (int i = 0; i < ri_nreads(*readinfo); ++i)
            hfile_write(outfile, &ri_readlen(*readinfo, i), uint64_t, 1);
        hfile_close_stream(outfile);
        anonymous_barrier;

        /// map-1
        clist = mapper(conf, ctx, readinfo, &runconf);
        anonymous_barrier;

        if (runconf.NNODES > 1) {
            lists = u128L1_partition(clist, &runconf);
            clist = shuffle(ctx, lists);
        }
        lists = u128L2_partition(clist, &runconf);
        void *idx = idx_create(lists, readinfo, ctx, &runconf);
        anonymous_barrier;

        align(conf, ctx, readinfo, idx, &runconf);
        anonymous_barrier;
    } else if (streq(argv[1], "assmbl")) {
        /// read data from temp file into clist
        hfile_t infile;
        hfile_init_input_stream(infile, "%s/alninfo-%03d.bin", hconf_tempdir(conf), hctx_mynode(ctx));
        u128_v blocksizes;
        hfile_read(infile, &blocksizes.n, size_t, 1);
        blocksizes.a = hmalloc(u128_t, blocksizes.n);
        hfile_read(infile, blocksizes.a, u128_t, blocksizes.n);
        blocksizes.m = blocksizes.n;
        clist = clist_create();
        hfile_close_stream(infile);
        hfile_init_input_stream(infile, "%s/aln-%03d.bin", hconf_tempdir(conf), hctx_mynode(ctx));
        for (int i = 0; i < blocksizes.n; ++i) {
            u_char *data = hmalloc(u_char, blocksizes.a[i].x);
            hfile_read(infile, data, u_char, blocksizes.a[i].x);
            clist_push_compressed(clist, data, blocksizes.a[i].x, blocksizes.a[i].y);
        }
        kv_destroy(blocksizes);
        LOG_INFO("Node %d: Clist blocks = %zd, Compfactor = %.2f.",
                 hctx_mynode(ctx), clist_length(clist), clist_compfactor(clist, ma_hit_t));

        if (runconf.NNODES > 1) {
            lists = hitL1_partition(clist, &runconf);
            clist = shuffle(ctx, lists);
        }

        lists = hitL2_partition(clist, &runconf);

        size_t nr; uint64_t rlen;
        hfile_init_input_stream(infile, "%s/rlen-%03d.bin", hconf_tempdir(conf), hctx_mynode(ctx));
        hfile_read(infile, &nr, size_t, 1);
        const uint32_t nseq = nr;
        sd_seq_t *seq = hmalloc(sd_seq_t, nseq);
        for (uint32_t i = 0; i < nseq; ++i) {
            hfile_read(infile, &rlen, uint64_t, 1);
            seq[i].id = i, seq[i].len = rlen;
        }
        hfile_close_stream(infile);

        /// select reads
        ma_sub_t *sub;
        uint32_t nseq_new = readselector(lists, seq, nseq, ctx, &runconf, &sub);
        anonymous_barrier;

        /// build graph
        build_graph(ctx, seq, nseq_new, sub, lists, &runconf, hconf_output(conf));
    } else LOG_ERR("Invalid phase %s", argv[1]);

    hctx_destroy(&ctx);
    hconf_destroy(&conf);
    return 0;
}