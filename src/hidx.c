//
// Created by sayan on 6/18/20.
//

#include "hiper3ga.h"
#include "kthread.h"
#include "hidx.h"

#include "khash.h"
#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(mm, uint32_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(mm) mmhash_t;
static int NNODES;
#define idx_key(k) (((k)/NNODES) >> PBITS)
#define idx_ptn(k) (((k)/NNODES) & PTNMASK)
static const char* HIDX_TAGSTR = "IDX4";
static int NIOTHREADS;
static bool SYNCREQ;
typedef struct {
    u32_v keys;
    u16_v counts;
    u64_v values;
    uint32_t *ptnsizes, *ptnoffsets;
    u128_t *mm; uint32_t nmm;
    u32_v r_offsets, r_counts;
    void *idx;
} request_t;
typedef struct {
    void *my_segment;
    void **remoteseg;           ///< array of segments dedicated for me by my peers
    size_t max_remote_segsize;
} channel_t;
typedef struct {
    uint32_t tag;
    int I, N, T;
    hctx_t *ctx;
    mmhash_t *h[NPTN];        ///< array of NPTN hash tables
    u64_v b[NPTN];            ///< position vector for minimizers appearing >1 times for each partition
    uint16_t max_occ;
    uintptr_t *peeridx;
    request_t *requests;
    channel_t *channels;
    runconf_t *runconf;
} idx_t;
static bool idx_is(void *self) {
    assert(self);
    return ((idx_t*)self)->tag == *((uint32_t *)HIDX_TAGSTR);
}
static inline bool channel_init(channel_t *c, int c_id, idx_t *idx) {
    c->remoteseg = hmalloc(void*, NNODES);
    size_t segsz = hctx_segsize(idx->ctx, 0);
    size_t max_segsize_per_node = SYNCREQ ? (segsz / 2) : (segsz / NNODES);
    c->max_remote_segsize = max_segsize_per_node / NIOTHREADS;
    for (int i = 0; i < NNODES; ++i) {
        void *segment = hctx_segaddr(idx->ctx, i);
        segment = SYNCREQ ? (segment + max_segsize_per_node) : (segment + max_segsize_per_node * idx->I);
        c->remoteseg[i] = segment + c->max_remote_segsize * c_id;
    }
    c->my_segment = SYNCREQ ? (hctx_segaddr(idx->ctx, idx->I) + c->max_remote_segsize * c_id) : c->remoteseg[idx->I];
}
static inline bool channel_clear(channel_t *c) {
    hfree(c->remoteseg);
}
static inline void request_init(request_t *r, idx_t *idx) {
    kv_init(r->keys); kv_init(r->counts), kv_init(r->values);
    kv_init(r->r_offsets); kv_init(r->r_counts);
    r->ptnsizes = hmalloc(uint32_t, NNODES);
    r->ptnoffsets = hmalloc(uint32_t, NNODES);
    r->idx = idx;
}
static inline void request_clear(request_t *r) {
    kv_destroy(r->keys); kv_destroy(r->counts); kv_destroy(r->values);
    kv_destroy(r->r_offsets); kv_destroy(r->r_counts);
    hfree(r->ptnsizes); hfree(r->ptnoffsets);
}

typedef struct {
    int ptnoff;
    clist_t **lists;
    idx_t *idx;
    u128_v *tbuf;
    readinfo_t *readinfo;
} shared_data_t;
static void indexer_fn(void *s, long i, int tid);
static const uint64_t *idx_get_local(idx_t *idx, uint32_t minier, uint16_t *n);
static void mm_idx_cal_max_occ_global(idx_t *idx, float f);
static void idx_get_from(int node, uint32_t *keys, size_t n_keys, uint16_t *counts, void *my_segment,
                         u64_v *v_values, idx_t *peeridx, void *peer_segment, size_t peer_segsize);

static idx_t *idx_setup(hctx_t *ctx, runconf_t *runconf) {
    LOG_INFO("Creating index..");
    /// setup
    idx_t *self = hmalloc(idx_t, 1);
    self->tag = *((uint32_t *)HIDX_TAGSTR);
    NNODES = self->N = hctx_num_nodes(ctx);
    self->T = runconf->NTHREADS;
    self->I = hctx_mynode(ctx);
    self->runconf = runconf;
    NIOTHREADS = hconf_niothreads(runconf->conf);
    if (NIOTHREADS < 0 || NIOTHREADS > self->T) NIOTHREADS = self->T;
    self->ctx = ctx;
    SYNCREQ = hconf_sync_requests();
    return self;
}
static void idx_sync(idx_t *self, hctx_t *ctx) {
    mm_idx_cal_max_occ_global(self, 0.001f);
    LOG_INFO("Max occ = %d.", self->max_occ);

    ulong * const my_segment = hctx_segaddr(ctx, -1);
    my_segment[0] = self->max_occ;
    my_segment[1] = (uintptr_t)self;
    anonymous_barrier;
    if (!self->I) { /// master
        for (uint i = 1; i < self->N; ++i)
            gasnet_get_nbi(my_segment + 2*i, i, hctx_segaddr(ctx, i), 2 * sizeof(ulong));
        gasnet_wait_syncnbi_gets();
        for (int i = 1; i < self->N; ++i)
            my_segment[0] = MAX(my_segment[0], my_segment[2*i]);
        for (int i = 1; i < self->N; ++i)
            my_segment[i+1] = my_segment[2*i+1];
        anonymous_barrier;
    } else {
        anonymous_barrier;
        gasnet_get(my_segment, 0, hctx_segaddr(ctx, 0), (self->N + 1) * sizeof(ulong));
    }
    self->peeridx = hmalloc(uintptr_t, self->N);
    self->max_occ = MIN(2048, my_segment[0]);
    for (int i = 0; i < self->N; ++i)
        self->peeridx[i] = my_segment[i+1];
    assert(self->peeridx[self->I] == (uintptr_t)self);
    if (!self->I) LOG_INFO("Global Max occ = %d.", self->max_occ);

    self->requests = hmalloc(request_t, self->T);
    for (int i = 0; i < self->T; ++i)
        request_init(&self->requests[i], self);

    self->channels = hmalloc(channel_t, NIOTHREADS);
    for (int i = 0; i < NIOTHREADS; ++i)
        channel_init(&self->channels[i], i, self);

    anonymous_barrier;
}

//////////// public functions
void *idx_create(clist_t** lists, readinfo_t *readinfo, hctx_t *ctx, runconf_t *runconf) {
    idx_t *self = idx_setup(ctx, runconf);

    u128_v tbuf[self->T];
    for (int i = 0; i < self->T; ++i) kv_init(tbuf[i]);
    shared_data_t sd = {.lists = lists, .idx = self, .tbuf = tbuf, .readinfo = readinfo, .ptnoff = 0};
    kt_forpool(runconf->tpool, indexer_fn, &sd, NPTN);
    for (int i = 0; i < self->T; ++i) kv_destroy(tbuf[i]);

    idx_sync(self, ctx);
    return self;
}
void idx_destroy(void **self_p) {
    assert(self_p && *self_p);
    assert(idx_is(*self_p));
    idx_t *self = *self_p;
    for (int i = 0; i < NPTN; ++i) {
        if (self->b[i].a) hfree(self->b[i].a);
        if (self->h[i]) kh_destroy(mm, self->h[i]);
    }
    for (int i = 0; i < self->T; ++i)
        request_clear(&self->requests[i]);
    hfree(self->requests);
    for (int i = 0; i < NIOTHREADS; ++i)
        channel_clear(&self->channels[i]);
    hfree(self->channels);
    hfree(self);
    *self_p = NULL;
}
#include "ksort.h"
#define sort_key_128L1(a) (((a).x)%NNODES)
KRADIX_SORT_INIT(128L1, u128_t, sort_key_128L1, 1)
void idx_setup_request(void *p_idx, int client_id, u128_t *mm, size_t nmm) {
    assert(idx_is(p_idx));
    idx_t *idx = (idx_t *)p_idx;

    /// sort by keys and then by partitions
    radix_sort_128L1(mm, mm + nmm);
    request_t *r = &idx->requests[client_id];

    memset(r->ptnsizes, 0, 4*NNODES);
    for (int i = 0; i < nmm; ++i)
        r->ptnsizes[sort_key_128L1(mm[i])]++;
    r->ptnoffsets[0] = 0;
    for (int i = 1; i < NNODES; ++i)
        r->ptnoffsets[i] = r->ptnoffsets[i-1] + r->ptnsizes[i-1];
    kv_reserve(uint32_t, r->keys, nmm);
    kv_reserve(uint16_t, r->counts, nmm);
    for (int i = 0; i < nmm; ++i) r->keys.a[i] = mm[i].x;
    r->values.n = 0;
    r->mm = mm, r->nmm = nmm;
#ifdef SANITY_CHECKS
    for (int j = r->ptnoffsets[idx->I]; j < r->ptnoffsets[idx->I] + r->ptnsizes[idx->I]; ++j)
        assert(r->keys.a[j] % NNODES == idx->I);
    for (int i = 0, k = 0; i < NNODES; ++i) {
        uint32_t *keys = r->keys.a + r->ptnoffsets[i];
        uint32_t n_keys = r->ptnsizes[i];
        for (int j = 0; j < n_keys; ++j, ++k) {
            assert(keys[j] % NNODES == i);
            assert(keys[j] == r->mm[k].x);
        }
    }
#endif
}
#include "utils.h"
typedef struct {
    idx_t *idx; int peer;
} shared_request_data_t;
static void sync_requester_fn(void *s, long c_id, int t_id) {
    shared_request_data_t *sd = (shared_request_data_t*)s;
    idx_t *idx = sd->idx;
    channel_t *c = &idx->channels[t_id];
    request_t *r = &idx->requests[c_id];
    int peer = sd->peer;
    idx_t *peeridx = (idx_t *) (idx->peeridx[peer]);
#ifdef SANITY_CHECKS
    for (int j = r->ptnoffsets[idx->I]; j < r->ptnoffsets[idx->I] + r->ptnsizes[idx->I]; ++j)
        assert(r->keys.a[j] % NNODES == idx->I);
    for (int i = 0, k = 0; i < NNODES; ++i) {
        uint32_t *keys = r->keys.a + r->ptnoffsets[i];
        uint32_t n_keys = r->ptnsizes[i];
        for (int j = 0; j < n_keys; ++j, ++k) {
            assert(keys[j] % NNODES == i);
            if(keys[j] != r->mm[k].x)
                LOG_ERR("Node %d, Thread %d: For client %ld at i=%d and j=%d, expected key %u, got %lu",
                        idx->I, t_id, c_id, i, j, keys[j], r->mm[k].x); /*r->mm[k].x = keys[j];*/
        }
    }
#endif
    if (r->ptnsizes[peer]) {
        idx_get_from(peer, r->keys.a + r->ptnoffsets[peer], r->ptnsizes[peer],
                     r->counts.a + r->ptnoffsets[peer], c->my_segment, &r->values,
                     peeridx, c->remoteseg[peer], c->max_remote_segsize);
    }
#ifdef SANITY_CHECKS
    for (int j = r->ptnoffsets[idx->I]; j < r->ptnoffsets[idx->I] + r->ptnsizes[idx->I]; ++j)
        assert(r->keys.a[j] % NNODES == idx->I);
    for (int i = 0, k = 0; i < NNODES; ++i) {
        uint32_t *keys = r->keys.a + r->ptnoffsets[i];
        uint32_t n_keys = r->ptnsizes[i];
        for (int j = 0; j < n_keys; ++j, ++k) {
            assert(keys[j] % NNODES == i);
            assert(keys[j] == r->mm[k].x);
        }
    }
#endif
}
static void requester_fn(void *p_idx, long client_id, int tid) {
    idx_t *idx = (idx_t *)p_idx;
    channel_t *c = &idx->channels[tid];
    request_t *r = &idx->requests[client_id];
#ifdef SANITY_CHECKS
    for (int j = r->ptnoffsets[idx->I]; j < r->ptnoffsets[idx->I] + r->ptnsizes[idx->I]; ++j)
        assert(r->keys.a[j] % NNODES == idx->I);
    for (int i = 0, k = 0; i < NNODES; ++i) {
        uint32_t *keys = r->keys.a + r->ptnoffsets[i];
        uint32_t n_keys = r->ptnsizes[i];
        for (int j = 0; j < n_keys; ++j, ++k) {
            assert(keys[j] % NNODES == i);
            assert(keys[j] == r->mm[k].x);
        }
    }
#endif
    for (int p = 1; p < NNODES; ++p) {
        int peer = (p + idx->I) % NNODES;
        idx_t *peeridx = (idx_t *) (idx->peeridx[peer]);
        if (!r->ptnsizes[peer]) continue;
        else
            idx_get_from(peer, r->keys.a + r->ptnoffsets[peer], r->ptnsizes[peer],
                         r->counts.a + r->ptnoffsets[peer], c->my_segment, &r->values,
                         peeridx, c->remoteseg[peer], c->max_remote_segsize);
    }
#ifdef SANITY_CHECKS
    for (int j = r->ptnoffsets[idx->I]; j < r->ptnoffsets[idx->I] + r->ptnsizes[idx->I]; ++j)
        assert(r->keys.a[j] % NNODES == idx->I);
    for (int i = 0, k = 0; i < NNODES; ++i) {
        uint32_t *keys = r->keys.a + r->ptnoffsets[i];
        uint32_t n_keys = r->ptnsizes[i];
        for (int j = 0; j < n_keys; ++j, ++k) {
            assert(keys[j] % NNODES == i);
            assert(keys[j] == r->mm[k].x);
        }
    }
#endif
}
void idx_wait(void *p_idx) {
    if (SYNCREQ) {
        assert(idx_is(p_idx));
        idx_t *idx = (idx_t *)p_idx;
        const int R = nearest_power2(idx->N);
        for (int round = 1; round < R; ++round) {
            gasnet_barrier_notify(round, 0);
            GASNET_SAFE(gasnet_barrier_wait(round, 0));
        }
    }
}
void idx_send_requests(void *p_idx) {
    assert(idx_is(p_idx));
    idx_t *idx = (idx_t *)p_idx;
    gasnett_local_mb();
    if (!SYNCREQ) kt_for(NIOTHREADS, requester_fn, p_idx, idx->T);
    else {
        shared_request_data_t sd = {.idx = idx};
        const int R = nearest_power2(idx->N);
        for (int round = 1; round < R; ++round) {
            int peer = idx->I ^round;
            gasnet_barrier_notify(round, 0);
            GASNET_SAFE(gasnet_barrier_wait(round, 0));
            if (peer < idx->N) {
                sd.peer = peer;
                kt_for(NIOTHREADS, sync_requester_fn, &sd, idx->T);
            }
        }
    }
    gasnett_local_mb();
}
#define REMOTE (((uint64_t)1)<<63u)
#define sort_key_addr(a) ((a).y)
KRADIX_SORT_INIT(mm_addr, u128_t, sort_key_addr, 8);
void idx_setup_response(void *p_idx, int client_id) {
    assert(idx_is(p_idx));
    idx_t *idx = (idx_t *)p_idx;
    request_t *r = &idx->requests[client_id];
#ifdef SANITY_CHECKS
    for (int i = 0, k = 0; i < NNODES; ++i) {
        uint32_t *keys = r->keys.a + r->ptnoffsets[i];
        uint32_t n_keys = r->ptnsizes[i];
        for (int j = 0; j < n_keys; ++j, ++k) {
            assert(keys[j] % NNODES == i);
            assert(keys[j] == r->mm[k].x);
        }
    }
#endif
    ulong offset = 0;
    if (SYNCREQ) {
        const int R = nearest_power2(idx->N);
        for (int round = 1; round < R; ++round) {
            int peer = idx->I ^round;
            if (peer < idx->N) {
                u128_t *a = r->mm + r->ptnoffsets[peer];
                size_t n = r->ptnsizes[peer];
                uint16_t *c = r->counts.a + r->ptnoffsets[peer];
                for (int j = 0; j < n; offset += c[j], ++j)
                    a[j].x = REMOTE | (offset << 16u) | c[j];
            }
        }
    } else {
        for (int i = 1; i < NNODES; ++i) {
            int peer = (i + idx->I) % NNODES;
            u128_t *a = r->mm + r->ptnoffsets[peer];
            size_t n = r->ptnsizes[peer];
            uint16_t *c = r->counts.a + r->ptnoffsets[peer];
            for (int j = 0; j < n; offset += c[j], ++j) {
                assert(sort_key_128L1(a[j]) == peer);
                assert(a[j].x == r->keys.a[r->ptnoffsets[peer] + j]);
                a[j].x = REMOTE | (offset << 16u) | c[j];
            }
        }
#ifdef SANITY_CHECKS
        for (int j = r->ptnoffsets[idx->I]; j < r->ptnoffsets[idx->I] + r->ptnsizes[idx->I]; ++j) {
            assert(r->mm[j].x == r->keys.a[j]);
            assert((r->mm[j].x & REMOTE) == 0);
        }
#endif
    }
    assert(offset < UINT32_MAX);

    /// sort by addresses, i.e., (read-id, position) tuples
    radix_sort_mm_addr(r->mm, r->mm + r->nmm);

    /// run-length encoding of read-ids
    r->r_offsets.n = r->r_counts.n = 0;
    offset = 0;
    unsigned i, j, count = 1;
    for (i = 0, j = 1; j < r->nmm; ++j) {
        if (mm_rid(r->mm[j]) == mm_rid(r->mm[i])) count++;
        else {
            kv_push(uint32_t, r->r_counts, count);
            kv_push(uint32_t, r->r_offsets, offset);
            i = j, offset += count, count = 1;
        }
    }
    kv_push(uint32_t, r->r_counts, count);
    kv_push(uint32_t, r->r_offsets, offset);
    assert(r->r_counts.n == r->r_offsets.n);
#ifdef SANITY_CHECKS
    for (i = 0; i < r->r_counts.n; ++i) {
        uint off = r->r_offsets.a[i], cnt = r->r_counts.a[i],
        rid = mm_rid(r->mm[off]), lastpos = mm_pos(r->mm[off]);
        u128_t *mm = r->mm + off; uint16_t ignore;
        for (j = 0; j < cnt; ++j) {
            assert(mm_rid(mm[j]) == rid);
            assert(mm_pos(mm[j]) >= lastpos);
            lastpos = mm_pos(mm[j]);
            if (!(mm[j].x & REMOTE)) {
                if (sort_key_128L1(mm[j]) != idx->I)
                    LOG_ERR("Node %d: Expected key %lu to be in node %lu.", idx->I, mm[j].x, sort_key_128L1(mm[j]));
                assert(idx_get_local(p_idx, mm[j].x, &ignore));
            }
        }
    }
#endif
}
unsigned idx_nreads(void *p_idx, int client_id) {
    assert(idx_is(p_idx));
    idx_t *idx = (idx_t *)p_idx;
    return idx->requests[client_id].r_counts.n;
}
u128_t * idx_get_mm(void *p_idx, int client_id, int i, uint *p_nmm, uint *p_rid) {
    assert(idx_is(p_idx));
    idx_t *idx = (idx_t *)p_idx;
    assert(i < idx->requests[client_id].r_counts.n);
    u128_t *mm = idx->requests[client_id].mm;
    uint offset = idx->requests[client_id].r_offsets.a[i];
    *p_nmm = idx->requests[client_id].r_counts.a[i];
    *p_rid = mm_rid(mm[offset]);
    return mm + offset;
}
const uint64_t * idx_get_vals(void *p_idx, int client_id, uint64_t key, uint16_t *n_vals) {
    if (key & REMOTE) {
        *n_vals = (uint16_t)key;
        uint32_t offset = key>>16;
        return ((idx_t*)p_idx)->requests[client_id].values.a + offset;
    } else return idx_get_local(p_idx, key, n_vals);
}

//////////// internal functions
#define u64_lt(a, b) ((a) < (b))
UTILS_INIT(64, uint64_t , u64_lt)
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, u128_t, sort_key_128x, 4)
static mmhash_t * build_hash(clist_t *list, size_t count, u128_t *mma, readinfo_t *readinfo, uint *p_bn, uint64_t **p_bp) {
    mmhash_t *h = NULL;
    uint bn = 0; ///size of the _p_ array in idx bucket
    uint64_t *bp = NULL; /// position array for minimizers appearing >1 times
    if (count) {
        int n, n_keys, max_n = 0;
        size_t j, start_a, start_p;

        /// decompress and sort
        long offset = 0, rc;
        while ((rc = clist_pop(list, u128_t, mma + offset, count - offset)) > 0)
            offset += rc;
        assert(offset == count);
        radix_sort_128x(mma, mma + count);

        /// transform read offsets into (read-id, position) pairs
        for (j = 0; j < count; ++j) {
            uint64_t rid = lower_bound_64(readinfo->cumulen, ri_nreads(*readinfo), mma[j].y >> 1) - 1;
            uint64_t pos = (mma[j].y >> 1) - readinfo->cumulen[rid];
            mma[j].y = rid << 32 | pos << 1 | (mma[j].y & 1);
        }

        /// count and preallocate
        for (j = 1, n = 1, n_keys = 0, bn = 0; j <= count; ++j) {
            if (j == count || mma[j].x != mma[j - 1].x) {
                ++n_keys;
                if (n > 1) bn += n;
                n = 1;
            } else ++n;
        }
        assert(n_keys);
        h = kh_init(mm);
        kh_resize(mm, h, n_keys);
        bp = hmalloc(uint64_t, bn);

        /// create the hash table
        for (j = 1, n = 1, start_a = start_p = 0; j <= count; ++j) {
            if (j == count || mma[j].x != mma[j - 1].x) {
                khint_t itr;
                int absent;
                u128_t *p = &mma[j - 1];
                itr = kh_put(mm, h, idx_key(p->x) << 1, &absent);
                assert(absent && j == start_a + n);
                if (n == 1) {
                    kh_key(h, itr) |= 1;
                    kh_val(h, itr) = p->y;
                } else {
                    int k;
                    for (k = 0; k < n; ++k)
                        bp[start_p + k] = mma[start_a + k].y;
                    kh_val(h, itr) = (uint64_t) start_p << 32 | n;
                    start_p += n;
                    if (n > max_n) max_n = n;
                }
                assert(kh_get(mm, h, idx_key(p->x) << 1) != kh_end(h));
                start_a = j, n = 1;
            } else ++n;
        }
    }

    *p_bn = bn, *p_bp = bp;
    return h;
}
static void indexer_fn(void *s, long i, int tid) {
    shared_data_t *sd = (shared_data_t*)s;
    idx_t *idx = sd->idx;
    clist_t *list = sd->lists[i];
    size_t count = clist_count(list);
    kv_reserve(u128_t, sd->tbuf[tid], count);
    u128_t *mma = sd->tbuf[tid].a;

    uint bn; ///size of the _p_ array in idx bucket
    uint64_t *bp; /// position array for minimizers appearing >1 times
    mmhash_t *h = build_hash(list, count, mma, sd->readinfo, &bn, &bp);
    clist_destroy(list);

    idx->b[i].n = idx->b[i].m = bn;
    idx->b[i].a = bp;
    idx->h[i] = h;
}
KSORT_INIT_GENERIC(uint32_t)
static void mm_idx_cal_max_occ_global(idx_t *idx, float f) {
    assert(idx_is(idx));
    int i;
    size_t n = 0;
    uint32_t thres;
    khint_t *a, k;
    uint32_t max_occ = UINT32_MAX;
    if (f <= 0.) return;
    for (i = 0; i < NPTN; ++i)
        if (idx->h[i]) n += kh_size(idx->h[i]);
    a = (uint32_t*)malloc(n * 4);
    for (i = n = 0; i < NPTN; ++i) {
        if (idx->h[i] == 0) continue;
        for (k = 0; k < kh_end(idx->h[i]); ++k) {
            if (!kh_exist(idx->h[i], k)) continue;
            a[n++] = kh_key(idx->h[i], k)&1? 1 : (uint32_t)kh_val(idx->h[i], k);
        }
    }
    thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
    free(a);
    idx->max_occ = (thres>UINT16_MAX) ? UINT16_MAX : thres;
}
static const uint64_t *idx_get_local(idx_t *idx, uint32_t minier, uint16_t *n) {
    assert(idx_is(idx));
    assert(minier % idx->N == idx->I);
    /// get the value for specified key from local hash tables
    const size_t p = idx_ptn(minier);
    const mmhash_t *h = idx->h[p];
    if (!h) LOG_ERR("This shouldn't be possible.");
    const uint64_t *bp = idx->b[p].a;
    if (!bp) LOG_ERR("This shouldn't be possible.");
    khint_t k = kh_get(mm, h, idx_key(minier)<<1);
    if (k == kh_end(h)) LOG_ERR("This shouldn't be possible. (key = %u)", minier);
    else if (kh_key(h, k)&1u) {
        *n = 1;
        return &kh_val(h, k);
    } else {
        *n = (uint32_t)kh_val(h, k);
        if (*n > idx->max_occ) *n = 0;
        return &bp[kh_val(h, k)>>32];
    }
}
static void idx_get_from(int node, uint32_t *keys, size_t n_keys, uint16_t *counts, void *my_segment,
                         u64_v *v_values, idx_t *peeridx, void *peer_segment, size_t peer_segsize) {
#ifdef SANITY_CHECKS
    for (int i = 0; i < n_keys; ++i)  assert(keys[i] % NNODES == node);
#endif
    if (n_keys * 12 > peer_segsize)
        LOG_ERR("Remote segment is too small (%zd MiB) for %zd keys.", peer_segsize<<20u, n_keys);

    gasnett_atomic_t flag = gasnett_atomic_init(0);
    GASNET_CALL(gasnet_AMRequestLong7, node, REQUEST_HANDLER_IDX, keys, n_keys * 4, peer_segment,
                SEND_PTR(peeridx), SEND_PTR(&flag), SEND_PTR(my_segment), (peer_segsize - n_keys*4)/8);
    unsigned n_vals;
    GASNET_BLOCKUNTIL((n_vals = gasnett_atomic_swap(&flag, 0, GASNETT_ATOMIC_ACQ)) > 0);
    gasnett_local_mb();
    uint16_t *my_countbuf = my_segment;
    uint64_t *my_valbuf = (uint64_t*)(((u_char*)my_segment) + n_keys * 4);
#ifdef SANITY_CHECKS
    size_t n_vals2 = 0;
    for (int x = 0; x < n_keys; ++x) n_vals2 += my_countbuf[x];
    assert(n_vals == n_vals2);
#endif
    memcpy(counts, my_countbuf, 2 * n_keys);
    kv_reserve(uint64_t, *v_values, v_values->n+n_vals);
    memcpy(v_values->a + v_values->n, my_valbuf, 8 * n_vals);
    v_values->n += n_vals;
}

////////////// AM handlers
void idx_request_handler(gasnet_token_t token, void *buf, size_t nbytes, PtrArg(idx), PtrArg(flag),
                         PtrArg(remotebuf), gasnet_handlerarg_t max_vals) {
    idx_t *idx = Ptr(idx);
    uint32_t *keys = buf;
    uint16_t *counts = buf;
    uint64_t *values = (uint64_t*)(((u_char*)buf) + nbytes);

    unsigned n_keys = nbytes/4, offset = 0, n_vals = 0;
    for (int i = 0; i < n_keys; ++i) {
        const uint64_t *v = idx_get_local(idx, keys[i], &counts[i]);
        n_vals += counts[i];
        if (!counts[i]) continue;
        if (offset + counts[i] > max_vals)
            LOG_ERR("# values (%u) exceeded buffer capacity (%u).", offset + counts[i], max_vals);
        else memcpy(values + offset, v, counts[i]*8);
        offset += counts[i];
    }
    assert(offset == n_vals);
    void *dest = Ptr(remotebuf);
    GASNET_CALL(gasnet_AMReplyLong3, token, RESPONSE_HANDLER_IDX,
                buf, nbytes + offset * 8, dest, PtrParam(flag), offset);
}

void idx_response_handler(gasnet_token_t token, void *buf, size_t nbytes, PtrArg(flag), gasnet_handlerarg_t n_vals) {
    gasnett_atomic_t *flag = Ptr(flag);
    gasnett_atomic_set(flag, n_vals, GASNETT_ATOMIC_REL);
}