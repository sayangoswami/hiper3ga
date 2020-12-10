//
// Created by sayan on 6/12/20.
//

#ifndef HIPER3GA_SHUFFLER_C
#define HIPER3GA_SHUFFLER_C

#include "hiper3ga.h"
#include "clist.h"
#include "utils.h"

clist_t * shuffle(hctx_t *ctx, clist_t **lists) {
    const int N = hctx_num_nodes(ctx);
    const int I = hctx_mynode(ctx);
    gasnett_atomic_t srcflag = gasnett_atomic_init(0);
    void *myseg = hctx_segaddr(ctx, I);
    clist_t *cl_out_p_copy = myseg;
    memcpy(cl_out_p_copy, lists[I], sizeof(clist_t));
    const size_t BATCHSZ = hctx_segsize(ctx, I) - sizeof(clist_t);
    size_t sent_sz = 0;
    const size_t mdt_count = gasnet_AMMaxMedium()/16;
    size_t sizes[mdt_count*2];
    const int R = nearest_power2(N);
    LOG_INFO("Node %d: Starting shuffle.", I);
    for (int round = 1; round < R; ++round) {
        if (!I) LOG_INFO("Round %d.", round);
        int peer = I ^ round;
        gasnet_barrier_notify(round, 0);
        GASNET_SAFE(gasnet_barrier_wait(round, 0));
        if (peer < N) {
            void *peeraddr = hctx_segaddr(ctx, peer);
            clist_t *destlist_p = peeraddr;
            void *destbuf = peeraddr + sizeof(clist_t);
            clist_t *srclist = lists[peer];
            size_t dst_offset = 0, count, size, blockcount = 0;
            u_char *data;

            while (clist_pop_compressed(srclist, &data, &size, &count)) {
                assert(count); assert(size);
                if (dst_offset + size >= BATCHSZ || blockcount == mdt_count) {
                    assert(dst_offset);
                    GASNET_CALL(gasnet_AMRequestMedium6, peer, SYNC_SEND_HANDLER_IDX, sizes, blockcount*16,
                                SEND_PTR(destbuf), SEND_PTR(destlist_p), SEND_PTR(&srcflag));
                    GASNET_BLOCKUNTIL(gasnett_atomic_swap(&srcflag, 0, GASNETT_ATOMIC_NONE) == 1);
                    dst_offset = 0;
                    blockcount = 0;
                }
                gasnet_put_bulk(peer, destbuf + dst_offset, data, size);
                sizes[2*blockcount] = count, sizes[2*blockcount+1] = size;
                blockcount++;
                dst_offset += size;
                sent_sz += size;
                hfree(data);
            }
            if (dst_offset) {
                assert(blockcount);
                GASNET_CALL(gasnet_AMRequestMedium6, peer, SYNC_SEND_HANDLER_IDX, sizes, blockcount*16,
                            SEND_PTR(destbuf), SEND_PTR(destlist_p), SEND_PTR(&srcflag));
                GASNET_BLOCKUNTIL(gasnett_atomic_swap(&srcflag, 0, GASNETT_ATOMIC_NONE) == 1);
            }
            clist_destroy(srclist);
        }
    }
    anonymous_barrier;
    LOG_INFO("Node %d: Sent %.3f GB of data to peers.", I, sent_sz/(1000000000.0));
    memcpy(lists[I], cl_out_p_copy, sizeof(clist_t));
    return lists[I];
}

void sync_send_handler(gasnet_token_t token, void *buf, size_t nbytes,
                       PtrArg(dstbuffer), PtrArg(dstlist), PtrArg(srcflag_p)) {
    clist_t *list = Ptr(dstlist);
    size_t *sizes = (size_t*)buf;
    u_char *srcdata = Ptr(dstbuffer);
    size_t blockcount = nbytes/16;
    size_t offset = 0;
    for (int i = 0; i < blockcount; ++i) {
        size_t count = sizes[2*i], csize = sizes[2*i+1];
        u_char *data = hmalloc(u_char, csize);
        memcpy(data, srcdata + offset, csize);
        clist_push_compressed(list, data, csize, count);
        offset += csize;
    }
    GASNET_CALL(gasnet_AMReplyShort2, token, ACK_HANDLER_IDX, PtrParam(srcflag_p));
}

void ack_handler(gasnet_token_t token, PtrArg(flag_p)) {
    gasnett_atomic_t *flag = Ptr(flag_p);
    gasnett_atomic_val_t oldval = gasnett_atomic_swap(flag, 1, GASNETT_ATOMIC_NONE);
    assert(oldval == 0);
}

#endif //HIPER3GA_SHUFFLER_C
