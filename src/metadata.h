//
// Created by sayan on 6/4/20.
//

#ifndef HIPER3GA_METADATA_H
#define HIPER3GA_METADATA_H

typedef struct {
    u128_v *reads;
    unsigned max_rdlen;
    uint64_t *cumulen;
    kvec_t(uint64_t) blockv;
    int blkoffset, nblocks;
} readinfo_t;

static readinfo_t *readinfo_create(hconf_t *conf, hctx_t *ctx, u128_v *reads) {
    const int I = hctx_mynode(ctx);
    const int N = hctx_num_nodes(ctx);
    readinfo_t *self = hmalloc(readinfo_t, 1);
    self->reads = reads;
    self->cumulen = hmalloc(uint64_t , reads->n);
    self->cumulen[0] = 0;
    self->max_rdlen = reads->a[0].y;
    for (int i = 1; i < reads->n; ++i) {
        self->cumulen[i] = self->cumulen[i-1] + reads->a[i-1].y;
        if (reads->a[i].y > self->max_rdlen) self->max_rdlen = reads->a[i].y;
    }
    const size_t max_nbases_per_block = hconf_max_bases_per_batch(conf);
    kv_init(self->blockv);
    size_t blksz = 0;
    uint32_t starting_read = 0, batch_readcount = 0;
    for (int i = 0; i < reads->n; ++i) {
        if (blksz + reads->a[i].y > max_nbases_per_block) {
            assert(batch_readcount);
            uint64_t block = ((uint64_t)starting_read)<<32u | batch_readcount;
            kv_push(uint64_t, self->blockv, block);
            starting_read += batch_readcount;
            batch_readcount = 0;
            blksz = 0;
        }
        ++batch_readcount;
        blksz += reads->a[i].y;
    }
    if (blksz) {
        assert(batch_readcount);
        uint64_t block = ((uint64_t)starting_read)<<32u | batch_readcount;
        kv_push(uint64_t, self->blockv, block);
        starting_read += batch_readcount;
        batch_readcount = 0;
        blksz = 0;
    }
    self->nblocks = self->blockv.n/N;
    if (I == N-1) self->nblocks += (self->blockv.n%N > 0);
    self->blkoffset = I * self->blockv.n/N;
    //const int my_lst_blk = my_1st_blk + my_nblocks;
    return self;
}
#define ri_nreads(readinfo) ((readinfo).reads->n)
#define ri_nblocks(readinfo) ((readinfo).nblocks)
#define ri_blkoffset(readinfo) ((readinfo).blkoffset)
#define ri_block(readinfo, bid) (readinfo).blockv.a[(bid)]
#define ri_block_startid(readinfo, bid) ((uint32_t)((readinfo).blockv.a[(bid)]>>32u))
#define ri_block_nreads(readinfo, bid) ((uint32_t)((readinfo).blockv.a[i]))
#define ri_readoffset(readinfo, rid) ((readinfo).reads->a[(rid)].x)
#define ri_readlen(readinfo, rid) ((readinfo).reads->a[(rid)].y)

#endif //HIPER3GA_METADATA_H
