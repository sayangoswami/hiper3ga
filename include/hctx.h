//
// Created by sayan on 9/13/19.
//

#ifndef HIPER3GA_HCTX_H
#define HIPER3GA_HCTX_H

hctx_t *hctx_create(int argc, char **argv);

void *hctx_destroy(hctx_t **self_p);

bool hctx_is(void *self);

int hctx_num_nodes(hctx_t *self);

int hctx_mynode(hctx_t *ctx);

/**
 * Segment size of peer (or self if peer = -1)
 * @param self
 * @param peer
 * @return
 */
size_t hctx_segsize(hctx_t *self, int peer);

/**
 * Segment address of peer (or self if peer = -1)
 * @param self
 * @param peer
 * @return
 */
void* hctx_segaddr(hctx_t *self, int peer);

#define SYNC_SEND_HANDLER_IDX 201
#define ACK_HANDLER_IDX 202
#define REQUEST_HANDLER_IDX 215
#define RESPONSE_HANDLER_IDX 216

#endif //HIPER3GA_HCTX_H
