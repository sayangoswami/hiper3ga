//
// Created by sayan on 9/13/19.
//

#include "hclasses.h"
static const char *TAGSTR = "CTXT";

struct _hctx_t {
    uint32_t tag;
    gasnet_node_t node_id;
    gasnet_node_t num_nodes;
    size_t segsize;                    ///< size of segment used for long active messages
    gasnet_seginfo_t *seginfo_table;   ///< an array of all remote segment addresses and sizes
    void *local_seg_addr;              ///< address of my registered segment accessible to peers
    gasnet_hsl_t lock;
};

extern void sync_send_handler(gasnet_token_t token, void *buf, size_t nbytes,
                              PtrArg(dstbuffer), PtrArg(dstlist), PtrArg(srcflag_p));
extern void ack_handler(gasnet_token_t token, PtrArg(flag_p));
extern void idx_request_handler(gasnet_token_t token, void *buf, size_t nbytes, PtrArg(idx), PtrArg(flag),
                                PtrArg(remotebuf), gasnet_handlerarg_t max_vals);
extern void idx_response_handler(gasnet_token_t token, void *buf, size_t nbytes, PtrArg(flag), gasnet_handlerarg_t n_vals);

hctx_t *hctx_create(int argc, char **argv) {
    /// bootstrap GASNet and initialize parameters
    GASNET_SAFE(gasnet_init(&argc, &argv));
    gasnett_backtrace_init(argv[0]);
    gasnet_set_waitmode(GASNET_WAIT_BLOCK);
    hctx_t *self = hmalloc(hctx_t, 1);
    self->tag = *((uint32_t *)TAGSTR);
    self->node_id = gasnet_mynode();
    self->num_nodes = gasnet_nodes();
    self->segsize = gasnet_getMaxGlobalSegmentSize();
    assert(self->segsize);
    size_t localsegsz = gasnet_getMaxLocalSegmentSize();

    /// initialize primary network resources
    gasnet_handlerentry_t htable[] = {
            {SYNC_SEND_HANDLER_IDX,   &sync_send_handler},
            {ACK_HANDLER_IDX,         &ack_handler},
            {REQUEST_HANDLER_IDX, &idx_request_handler},
            {RESPONSE_HANDLER_IDX, &idx_response_handler}
    };
    GASNET_SAFE(gasnet_attach(htable, sizeof(htable) / sizeof(gasnet_handlerentry_t), self->segsize, 0));
    if (self->node_id == 0) {
        LOG_INFO("Max local segsize = %zd bytes.", localsegsz);
        LOG_INFO("Attached segsize = %zd bytes.", self->segsize);
        LOG_INFO("Max args = %zd, Max medium = %zd bytes, Max long request / reply = %zd / %zd bytes.",
                 gasnet_AMMaxArgs(), gasnet_AMMaxMedium(), gasnet_AMMaxLongRequest(), gasnet_AMMaxLongReply());
    }
    self->seginfo_table = hmalloc(gasnet_seginfo_t, self->num_nodes);
    GASNET_SAFE(gasnet_getSegmentInfo(self->seginfo_table, self->num_nodes));

    self->local_seg_addr = self->seginfo_table[self->node_id].addr;
    assert(self->segsize == self->seginfo_table[self->node_id].size);
    gasnet_hsl_init(&self->lock);
    anonymous_barrier;

    return self;
}

void *hctx_destroy(hctx_t **self_p) {
    if (self_p && *self_p) {
        hctx_t *self = *self_p;
        assert(hctx_is(self));
        hfree(self->seginfo_table);
        gasnet_hsl_destroy(&self->lock);
        hfree(self);
        *self_p = NULL;
        anonymous_barrier;
        gasnet_exit(0);
    }
}

bool hctx_is(void *self) {
    assert(self);
    return ((hctx_t*)self)->tag == *((uint32_t *)TAGSTR);
}

int hctx_num_nodes(hctx_t *self) {
    assert(hctx_is(self));
    return self->num_nodes;
}

int hctx_mynode(hctx_t *self) {
    assert(hctx_is(self));
    return self->node_id;
}

size_t hctx_segsize(hctx_t *self, int peer) {
    assert(hctx_is(self));
    if (peer >= 0) return self->seginfo_table[peer].size;
    else return self->segsize;
}

void* hctx_segaddr(hctx_t *self, int peer) {
    assert(hctx_is(self));
    if (peer >= 0) return self->seginfo_table[peer].addr;
    else return self->local_seg_addr;
}
