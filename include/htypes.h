//
// Created by sayan on 9/24/19.
//

#ifndef HIPER3GA_HTYPES_H
#define HIPER3GA_HTYPES_H

typedef struct {uint64_t x, y;} u128_t;
typedef kvec_t(u128_t) u128_v;
typedef kvec_t(uint64_t) u64_v;
typedef kvec_t(uint32_t) u32_v;
typedef kvec_t(uint16_t) u16_v;

#define u128_leq(a, b) (        \
    ((a).x < (b).x) ? 1 : (     \
        ((a).x > (b).x) ? 0 :   \
            ((a).y <= (b).y)    \
    )                           \
)
#define u128_lt(a, b) (         \
    ((a).x < (b).x) ? 1 : (     \
        ((a).x > (b).x) ? 0 :   \
            ((a).y < (b).y)     \
    )                           \
)
#define u128_geq(a, b) (        \
    ((a).x > (b).x) ? 1 : (     \
        ((a).x < (b).x) ? 0 :   \
            ((a).y >= (b).y)    \
    )                           \
)
#define u128_gt(a, b) (         \
    ((a).x > (b).x) ? 1 : (     \
        ((a).x < (b).x) ? 0 :   \
            ((a).y > (b).y)     \
    )                           \
)
#define u128_min(a, b) ((u128_lt((a), (b))) ? (a) : (b))
#define u128_cmp(a, b) ((u128_lt((a), (b))) ? -1 : (u128_gt((a), (b))))

#define range_t(type) struct { type *begin, *end; }
#define range_init(data_array, data_array_size, range_array, range_array_size) do {     \
    size_t na = data_array_size/range_array_size;                                       \
    range_array[0].begin = data_array;                                                  \
    range_array[range_array_size-1].end = data_array + data_array_size;                 \
    for (int i = 1; i < range_array_size; ++i) {                                        \
        range_array[i-1].end = range_array[i].begin = range_array[i-1].begin + na;      \
    }                                                                                   \
} while(0)
#define range_size(r) ((r).end - (r).begin)
typedef range_t(uint64_t) u64_r;

typedef struct {
    void *tpool;
    int NTHREADS;
    int NNODES;
    pthread_barrier_t *barrier;
    char log_prefix[64];
    hconf_t *conf;
} runconf_t;

//typedef struct {
//    uint32_t key;
//    uint32_t _data[2]; /// count (14 bits), block id (20), block offset (30 bits)
//} triple_t;
//#define t_0                   {0, {0,0}}
//#define t_key(t)              ((t).key)
//#define t_value(t)            (((uint64_t*)((t)._data))[0])
//#define v_count(v)            (((uint64_t)(v))>>50u)
//#define v_blkid(v)            ((((uint64_t)(v))<<14u)>>44u)
//#define v_blkoffset(v)        ((((uint64_t)(v))<<34u)>>34u)
//#define t_getcount(t)         v_count(t_value(t))
//#define t_set(t, k, c, i, o)  (t).key = (k), t_value(t) = ((uint64_t)(c))<<50u | ((uint64_t)(i))<<30u | (uint64_t)(o)

typedef struct {
    uint32_t cnt:31, rev:1;
    uint32_t rid:31, rep:1;
    uint32_t len;
    int32_t qs, qe, rs, re;
} mm_reg1_t;

typedef struct {
    uint32_t qid, ql, qs, qe, tid, tl, ts, te;
    uint32_t ml:31, rev:1, bl;
} paf_rec_t;

typedef kvec_t(paf_rec_t) paf_rec_v;

typedef struct {
    uint64_t qns;
    uint32_t qe, tn, ts, te;
    uint32_t ml:31, rev:1;
    uint32_t bl:31, del:1;
} ma_hit_t;
typedef kvec_t(ma_hit_t) ma_hit_v;

typedef struct {
    /*char *name;*/
    uint32_t id, len, aux:31, del:1;
} sd_seq_t;

typedef struct {
    uint32_t s:31, del:1, e;
} ma_sub_t;

typedef struct {
    uint64_t ul;
    uint32_t v;
    uint32_t ol:31, del:1;
} asg_arc_t;

#define MA_HT_INT        (-1)
#define MA_HT_QCONT      (-2)
#define MA_HT_TCONT      (-3)
#define MA_HT_SHORT_OVLP (-4)

static int ma_hit2arc(const ma_hit_t *h, int ql, int tl, int _max_hang, float _int_frac, int _min_ovlp, asg_arc_t *p)
{
    int32_t tl5, tl3, ext5, ext3, qs = (int32_t)h->qns;
    uint32_t u, v, l; // u: query end; v: target end; l: length from u to v
    if (h->rev) tl5 = tl - h->te, tl3 = h->ts; // tl5: 5'-end overhang (on the query strand); tl3: similar
    else tl5 = h->ts, tl3 = tl - h->te;
    ext5 = qs < tl5? qs : tl5;
    ext3 = ql - h->qe < tl3? ql - h->qe : tl3;
    if (ext5 > _max_hang || ext3 > _max_hang || h->qe - qs < (h->qe - qs + ext5 + ext3) * _int_frac)
        return MA_HT_INT;
    if (qs <= tl5 && ql - h->qe <= tl3) return MA_HT_QCONT; // query contained
    else if (qs >= tl5 && ql - h->qe >= tl3) return MA_HT_TCONT; // target contained
    else if (qs > tl5) u = 0, v = !!h->rev, l = qs - tl5;
    else u = 1, v = !h->rev, l = (ql - h->qe) - tl3;
    if (h->qe - qs + ext5 + ext3 < _min_ovlp || h->te - h->ts + ext5 + ext3 < _min_ovlp) return MA_HT_SHORT_OVLP; // short overlap
    u |= h->qns>>32<<1, v |= h->tn<<1;
    p->ul = (uint64_t)u<<32 | l, p->v = v, p->ol = ql - l, p->del = 0;
    return l;
}

#endif //HIPER3GA_HTYPES_H
