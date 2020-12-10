//
// Created by sayan on 6/8/20.
//

#include "hiper3ga.h"
#include "kthread.h"
#include "ksort.h"
#include "clist.h"

typedef struct {
    uint32_t len:31, del:1;
} asg_seq_t;

typedef struct {
    uint32_t m_arc, n_arc:31, is_srt:1;
    asg_arc_t *arc;
    uint32_t m_seq, n_seq:31, is_symm:1;
    asg_seq_t *seq;
    uint64_t *idx;
} asg_t;
typedef struct { size_t n, m; uint64_t *a; } asg64_v;

typedef struct {
    uint32_t len:31, circ:1; ///< len: length of the unitig; circ: circular if non-zero
    uint32_t start, end; /// start: starting vertex in the string graph; end: ending vertex
    uint32_t m, n; /// number of reads
    uint64_t *a; /// list of reads
    char *s; /// unitig sequence is not null
} ma_utg_t;
typedef kvec_t(ma_utg_t) ma_utg_v;

typedef struct {
    ma_utg_v u;
    asg_t *g;
} ma_ug_t;

#define asg_arc_len(arc) ((uint32_t)(arc).ul)
#define asg_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define asg_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])

#define asg_arc_key(a) ((a).ul)
KRADIX_SORT_INIT(asg, asg_arc_t, asg_arc_key, 8)

static int asg_arc_del_trans(asg_t *g, int fuzz);
static int asg_cut_tip(asg_t *g, int max_ext);
static int asg_pop_bubble(asg_t *g, int max_dist);
static int asg_cut_internal(asg_t *g, int max_ext);
static int asg_cut_biloop(asg_t *g, int max_ext);
static int asg_arc_del_short(asg_t *g, float drop_ratio);
static void asg_cleanup(asg_t *g);
static void asg_seq_set(asg_t *g, int sid, int len, int del);
static inline asg_arc_t *asg_arc_pushp(asg_t *g);
static ma_ug_t *ma_ug_gen(asg_t *g);
static void ma_ug_print(const ma_ug_t *ug, const sd_seq_t *seq, int nseq, const ma_sub_t *sub, FILE *fp);

typedef struct {
    const sd_seq_t *seq; /// lengths of reads (readonly)
    uint32_t n_seq; /// number of reads
    const ma_sub_t *sub; /// readonly
    clist_t **lists; /// array of list of hits
    hbitvec_t *bv; /// bit vector for thread-safe writes
    asg_arc_t *arcs[ASM_NPTN];
    size_t num_arcs[ASM_NPTN];
    ma_hit_v *hitbuf;
} shared_data_t;
static void builder_fn(void *s, long id, int tid) {
    shared_data_t *sd = (shared_data_t*)s;
    ma_hit_v *hitbuf = &sd->hitbuf[tid];
    hitbuf->n = clist_count(sd->lists[id]); /// number of hits
    kv_reserve(ma_hit_t, *hitbuf, hitbuf->n);   /// array of hits
    long offset = 0, count;
    while ((count = clist_pop(sd->lists[id], ma_hit_t, hitbuf->a + offset, hitbuf->n - offset)) != 0)
        offset += count;
    assert(offset == hitbuf->n);
    clist_destroy(sd->lists[id]);
    kvec_t(asg_arc_t) arcs = {0, 0, 0};
    size_t i;
    for (i = 0; i < hitbuf->n; ++i) {
        int r;
        asg_arc_t t; //, *p;
        const ma_hit_t *h = &hitbuf->a[i];
        uint32_t qn = h->qns>>32;
        int ql = sd->sub? sd->sub[qn].e - sd->sub[qn].s : sd->seq[qn].len;
        int tl = sd->sub? sd->sub[h->tn].e - sd->sub[h->tn].s : sd->seq[h->tn].len;
        r = ma_hit2arc(h, ql, tl, max_hang, int_frac, min_ovlp, &t);
        if (r >= 0) {
            if (qn == h->tn) { /// self match
                if ((uint32_t)h->qns == h->ts && h->qe == h->te && h->rev)
                    hbitvec_set(sd->bv, qn);
                continue;
            }
            kv_push(asg_arc_t, arcs, t);
        } else if (r == MA_HT_QCONT)
            hbitvec_set(sd->bv, qn);
    }
    sd->arcs[id] = arcs.a; sd->num_arcs[id] = arcs.n;
}

void build_graph(hctx_t *ctx, const sd_seq_t *seq, uint32_t n_seq,
                 const ma_sub_t *sub, clist_t ** lists, runconf_t *runconf, char *output) {

    const int I = hctx_mynode(ctx);
    const int N = hctx_num_nodes(ctx);
    ma_hit_v hitbuf[runconf->NTHREADS];
    for (int i = 0; i < runconf->NTHREADS; ++i) kv_init(hitbuf[i]);
    shared_data_t sd = {
            .seq = seq, .n_seq = n_seq, .sub = sub, .lists = lists, .bv = hbitvec_new(n_seq), .hitbuf = hitbuf
    };
    kt_forpool(runconf->tpool, builder_fn, &sd, ASM_NPTN);
    size_t total_num_arcs = 0;
    for (int i = 0; i < ASM_NPTN; ++i) total_num_arcs += sd.num_arcs[i];
    LOG_INFO("Node %d: %zd arcs created.", I, total_num_arcs);

    ///////////////// collate arcs
    size_t *narcs_p = hctx_segaddr(ctx, I);
    asg_arc_t *allarcs = hctx_segaddr(ctx, I) + sizeof(size_t);
    *narcs_p = total_num_arcs;
    assert(total_num_arcs * sizeof(asg_arc_t) + sizeof(size_t) < hctx_segsize(ctx, I));
    for (size_t i = 0, off = 0; i < ASM_NPTN; ++i) {
        memcpy(allarcs+off, sd.arcs[i], sd.num_arcs[i] * sizeof(asg_arc_t));
        hfree(sd.arcs[i]);
        off += sd.num_arcs[i];
    }
    anonymous_barrier;

    if (N > 1) {
        if (I == 0) {
            size_t rcounts[N], total_count = 0;
            gasnet_valget_handle_t cnthandles[N];

            for (uint i = 0; i < N; ++i) cnthandles[i] = gasnet_get_nb_val(i, hctx_segaddr(ctx, i), sizeof(size_t));
            for (uint i = 0; i < N; ++i) {
                rcounts[i] = gasnet_wait_syncnb_valget(cnthandles[i]);
                total_count += rcounts[i];
            } assert(total_count * sizeof(asg_arc_t) + sizeof(size_t) < hctx_segsize(ctx, 0));

            size_t off = rcounts[0];
            for (uint i = 1; i < N; ++i) {
                gasnet_get_nbi(allarcs + off, i, hctx_segaddr(ctx, i) + sizeof(size_t), rcounts[i] * sizeof(asg_arc_t));
                off += rcounts[i];
            } gasnet_wait_syncnbi_gets();
            assert(off == total_count);
            *narcs_p = total_count;
        }
        anonymous_barrier;
        size_t count = 0;
        uint32_t *array = hbitvec_raw(sd.bv, &count);
        assert(count == n_seq);
        const size_t array_size = count/32 + (count%32 != 0);
        if (I) memcpy(hctx_segaddr(ctx, I), array, array_size * sizeof(uint32_t));
        anonymous_barrier;
        if (I == 0) {
            uint32_t *remote_array = hmalloc(uint32_t, array_size);
            for (int n = 1; n < N; ++n) {
                gasnet_get_bulk(remote_array, n, hctx_segaddr(ctx, n), array_size * sizeof(uint32_t));
                for (size_t i = 0; i < array_size; ++i) array[i] |= remote_array[i];
            }
            hfree(remote_array);
        }
    }

    //////////////////////////////////////

    if (I == 0) {
        /// TODO - transitively reduce edges and compress
        asg_t *g = hmalloc(asg_t, 1);
        for (int i = 0; i < n_seq; ++i) {
            if (sub) asg_seq_set(g, i, sub[i].e - sub[i].s, (sub[i].del || seq[i].del));
            else asg_seq_set(g, i, seq[i].len, seq[i].del);
        }
        size_t narcs = *narcs_p;
        for (int i = 0; i < narcs; ++i) {
            asg_arc_t *p = asg_arc_pushp(g);
            *p = allarcs[i];
        }
        for (size_t i = 0; i < n_seq; ++i) if (hbitvec_get(sd.bv, i)) g->seq[i].del = 1;
        asg_cleanup(g);
        LOG_INFO("Node %d: Read %d arcs.", I, g->n_arc);

        LOG_INFO("[M::%s] ===> Step 4.1: transitive reduction <===", __func__);
        asg_arc_del_trans(g, gap_fuzz);

        LOG_INFO("[M::%s] ===> Step 4.2: initial tip cutting and bubble popping <===", __func__);
        asg_cut_tip(g, max_ext);
        asg_pop_bubble(g, bub_dist);

        LOG_INFO("[M::%s] ===> Step 4.3: cutting short overlaps (%d rounds in total) <===", __func__, n_rounds + 1);
        for (int i = 0; i <= n_rounds; ++i) {
            float r = min_ovlp_drop_ratio + (max_ovlp_drop_ratio - min_ovlp_drop_ratio) / n_rounds * i;
            if (asg_arc_del_short(g, r) != 0) {
                asg_cut_tip(g, max_ext);
                asg_pop_bubble(g, bub_dist);
            }
        }

        LOG_INFO("[M::%s] ===> Step 4.4: removing short internal sequences and bi-loops <===", __func__);
        asg_cut_internal(g, 1);
        asg_cut_biloop(g, max_ext);
        asg_cut_tip(g, max_ext);
        asg_pop_bubble(g, bub_dist);

        LOG_INFO("[M::%s] ===> Step 4.5: aggressively cutting short overlaps <===", __func__);
        if (asg_arc_del_short(g, final_ovlp_drop_ratio) != 0) {
            asg_cut_tip(g, max_ext);
            asg_pop_bubble(g, bub_dist);
        }

        LOG_INFO("[M::%s] ===> Step 5: generating unitigs <===", __func__);
        ma_ug_t *ug = ma_ug_gen(g);
        FILE *fp = fopen(output, "w");
        if(!fp) LOG_ERR("Could not open file %s because %s.", output, strerror(errno));
        ma_ug_print(ug, seq, n_seq, sub, fp);
        fclose(fp);

    } else {
        /// nothing to do
    }
    anonymous_barrier;
}

static void asg_seq_set(asg_t *g, int sid, int len, int del) {
    if (sid >= g->m_seq) {
        g->m_seq = sid + 1;
        kv_roundup32(g->m_seq);
        g->seq = hrealloc(g->seq, asg_seq_t, g->m_seq);
    }
    if (sid >= g->n_seq) g->n_seq = sid + 1;
    g->seq[sid].len = len;
    g->seq[sid].del = !!del;
}

static inline asg_arc_t *asg_arc_pushp(asg_t *g) {
    if (g->n_arc == g->m_arc) {
        g->m_arc = g->m_arc? g->m_arc<<1 : 16;
        g->arc = hrealloc(g->arc, asg_arc_t, g->m_arc);
    }
    return &g->arc[g->n_arc++];
}

/// hard remove arcs marked as "del"
static void asg_arc_rm(asg_t *g) {
    uint32_t e, n;
    for (e = n = 0; e < g->n_arc; ++e) {
        uint32_t u = g->arc[e].ul>>32, v = g->arc[e].v;
        if (!g->arc[e].del && !g->seq[u>>1].del && !g->seq[v>>1].del)
            g->arc[n++] = g->arc[e];
    }
    if (n < g->n_arc) { // arc index is out of sync
        if (g->idx) hfree(g->idx);
        g->idx = 0;
    }
    g->n_arc = n;
}

static void asg_arc_sort(asg_t *g) {
    radix_sort_asg(g->arc, g->arc + g->n_arc);
}

static uint64_t *asg_arc_index_core(size_t max_seq, size_t n, const asg_arc_t *a) {
    size_t i, last;
    uint64_t *idx;
    idx = hmalloc(uint64_t, max_seq*2);
    for (i = 1, last = 0; i <= n; ++i)
        if (i == n || a[i-1].ul>>32 != a[i].ul>>32)
            idx[a[i-1].ul>>32] = (uint64_t)last<<32 | (i - last), last = i;
    return idx;
}

static void asg_arc_index(asg_t *g) {
    if (g->idx) hfree(g->idx);
    g->idx = asg_arc_index_core(g->n_seq, g->n_arc, g->arc);
}

static void asg_cleanup(asg_t *g) {
    asg_arc_rm(g);
    if (!g->is_srt) {
        asg_arc_sort(g);
        g->is_srt = 1;
    }
    if (g->idx == 0) asg_arc_index(g);
}

/// delete multi-arcs
static int asg_arc_del_multi(asg_t *g)
{
    uint32_t *cnt, n_vtx = g->n_seq * 2, n_multi = 0, v;
    cnt = hmalloc(uint32_t, n_vtx);
    for (v = 0; v < n_vtx; ++v) {
        asg_arc_t *av = asg_arc_a(g, v);
        int32_t i, nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        for (i = nv - 1; i >= 0; --i) ++cnt[av[i].v];
        for (i = nv - 1; i >= 0; --i)
            if (--cnt[av[i].v] != 0)
                av[i].del = 1, ++n_multi;
    }
    hfree(cnt);
    if (n_multi) asg_cleanup(g);
    LOG_INFO("[M::%s] removed %d multi-arcs", __func__, n_multi);
    return n_multi;
}

/// remove asymmetric arcs: u->v is present, but v'->u' not
static int asg_arc_del_asymm(asg_t *g) {
    uint32_t e, n_asymm = 0;
    for (e = 0; e < g->n_arc; ++e) {
        uint32_t v = g->arc[e].v^1, u = g->arc[e].ul>>32^1;
        uint32_t i, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        for (i = 0; i < nv; ++i)
            if (av[i].v == u) break;
        if (i == nv) g->arc[e].del = 1, ++n_asymm;
    }
    if (n_asymm) asg_cleanup(g);
    LOG_INFO("[M::%s] removed %d asymmetric arcs", __func__, n_asymm);
    return n_asymm;
}

static void asg_symm(asg_t *g) {
    asg_arc_del_multi(g);
    asg_arc_del_asymm(g);
    g->is_symm = 1;
}

/// transitive reduction; see Myers, 2005
static int asg_arc_del_trans(asg_t *g, int fuzz) {
    uint8_t *mark;
    uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;

    mark = hmalloc(uint8_t, n_vtx);
    for (v = 0; v < n_vtx; ++v) {
        uint32_t L, i, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (nv == 0) continue; // no hits
        if (g->seq[v>>1].del) {
            for (i = 0; i < nv; ++i) av[i].del = 1, ++n_reduced;
            continue;
        }
        for (i = 0; i < nv; ++i) mark[av[i].v] = 1;
        L = asg_arc_len(av[nv-1]) + fuzz;
        for (i = 0; i < nv; ++i) {
            uint32_t w = av[i].v;
            uint32_t j, nw = asg_arc_n(g, w);
            asg_arc_t *aw = asg_arc_a(g, w);
            if (mark[av[i].v] != 1) continue;
            for (j = 0; j < nw && asg_arc_len(aw[j]) + asg_arc_len(av[i]) <= L; ++j)
                if (mark[aw[j].v]) mark[aw[j].v] = 2;
        }
#if 0
        for (i = 0; i < nv; ++i) {
                        uint32_t w = av[i].v;
                        uint32_t j, nw = asg_arc_n(g, w);
                        asg_arc_t *aw = asg_arc_a(g, w);
                        for (j = 0; j < nw && (j == 0 || asg_arc_len(aw[j]) < fuzz); ++j)
                                if (mark[aw[j].v]) mark[aw[j].v] = 2;
                }
#endif
        for (i = 0; i < nv; ++i) {
            if (mark[av[i].v] == 2) av[i].del = 1, ++n_reduced;
            mark[av[i].v] = 0;
        }
    }
    hfree(mark);
    LOG_INFO("[M::%s] transitively reduced %d arcs", __func__, n_reduced);
    if (n_reduced) {
        asg_cleanup(g);
        asg_symm(g);
    }
    return n_reduced;
}

/// set asg_arc_t::del for v->w
static inline void asg_arc_del(asg_t *g, uint32_t v, uint32_t w, int del) {
    uint32_t i, nv = asg_arc_n(g, v);
    asg_arc_t *av = asg_arc_a(g, v);
    for (i = 0; i < nv; ++i)
        if (av[i].v == w) av[i].del = !!del;
}

#define ASG_ET_MERGEABLE 0
#define ASG_ET_TIP       1
#define ASG_ET_MULTI_OUT 2
#define ASG_ET_MULTI_NEI 3

static inline int asg_is_utg_end(const asg_t *g, uint32_t v, uint64_t *lw) {
    uint32_t w, nv, nw, nw0, nv0 = asg_arc_n(g, v^1);
    int i, i0 = -1;
    asg_arc_t *aw, *av = asg_arc_a(g, v^1);
    for (i = nv = 0; i < nv0; ++i)
        if (!av[i].del) i0 = i, ++nv;
    if (nv == 0) return ASG_ET_TIP; // tip
    if (nv > 1) return ASG_ET_MULTI_OUT; // multiple outgoing arcs
    if (lw) *lw = av[i0].ul<<32 | av[i0].v;
    w = av[i0].v ^ 1;
    nw0 = asg_arc_n(g, w);
    aw = asg_arc_a(g, w);
    for (i = nw = 0; i < nw0; ++i)
        if (!aw[i].del) ++nw;
    if (nw != 1) return ASG_ET_MULTI_NEI;
    return ASG_ET_MERGEABLE;
}

static int asg_extend(const asg_t *g, uint32_t v, int max_ext, asg64_v *a) {
    int ret;
    uint64_t lw;
    a->n = 0;
    kv_push(uint64_t, *a, v);
    do {
        ret = asg_is_utg_end(g, v^1, &lw);
        if (ret != 0) break;
        kv_push(uint64_t, *a, lw);
        v = (uint32_t)lw;
    } while (--max_ext > 0);
    return ret;
}

/// set asg_arc_t::del and asg_seq_t::del to 1 for sequence s and all its associated arcs
static inline void asg_seq_del(asg_t *g, uint32_t s) {
    uint32_t k;
    g->seq[s].del = 1;
    for (k = 0; k < 2; ++k) {
        uint32_t i, v = s<<1 | k;
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        for (i = 0; i < nv; ++i) {
            av[i].del = 1;
            asg_arc_del(g, av[i].v^1, v^1, 1);
        }
    }
}

static int asg_cut_tip(asg_t *g, int max_ext) {
    asg64_v a = {0, 0, 0};
    uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v >> 1].del) continue;
        if (asg_is_utg_end(g, v, 0) != ASG_ET_TIP) continue; // not a tip
        if (asg_extend(g, v, max_ext, &a) == ASG_ET_MERGEABLE) continue; // not a short unitig
        for (i = 0; i < a.n; ++i)
            asg_seq_del(g, (uint32_t) a.a[i] >> 1);
        ++cnt;
    }
    hfree(a.a);
    if (cnt > 0) asg_cleanup(g);
    LOG_INFO("[M::%s] cut %d tips", __func__, cnt);
    return cnt;
}

typedef struct {
    uint32_t p; // the optimal parent vertex
    uint32_t d; // the shortest distance from the initial vertex
    uint32_t c; // max count of reads
    uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
} binfo_t;

typedef struct {
    binfo_t *a;
    kvec_t(uint32_t) S; // set of vertices without parents
    kvec_t(uint32_t) T; // set of tips
    kvec_t(uint32_t) b; // visited vertices
    kvec_t(uint32_t) e; // visited edges/arcs
} buf_t;

/// in a resolved bubble, mark unused vertices and arcs as "reduced"
static void asg_bub_backtrack(asg_t *g, uint32_t v0, buf_t *b) {
    uint32_t i, v;
    assert(b->S.n == 1);
    for (i = 0; i < b->b.n; ++i)
        g->seq[b->b.a[i]>>1].del = 1;
    for (i = 0; i < b->e.n; ++i) {
        asg_arc_t *a = &g->arc[b->e.a[i]];
        a->del = 1;
        asg_arc_del(g, a->v^1, a->ul>>32^1, 1);
    }
    v = b->S.a[0];
    do {
        uint32_t u = b->a[v].p; // u->v
        g->seq[v>>1].del = 0;
        asg_arc_del(g, u, v, 0);
        asg_arc_del(g, v^1, u^1, 0);
        v = u;
    } while (v != v0);
}

/// count the number of outgoing arcs, excluding reduced arcs
static inline int count_out(const asg_t *g, uint32_t v) {
    uint32_t i, n, nv = asg_arc_n(g, v);
    const asg_arc_t *av = asg_arc_a(g, v);
    for (i = n = 0; i < nv; ++i)
        if (!av[i].del) ++n;
    return n;
}

/// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
static uint64_t asg_bub_pop1(asg_t *g, uint32_t v0, int max_dist, buf_t *b) {
    uint32_t i, n_pending = 0;
    uint64_t n_pop = 0;
    if (g->seq[v0>>1].del) return 0; // already deleted
    if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
    b->S.n = b->T.n = b->b.n = b->e.n = 0;
    b->a[v0].c = b->a[v0].d = 0;
    kv_push(uint32_t, b->S, v0);
    do {
        uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
        uint32_t nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        assert(nv > 0);
        for (i = 0; i < nv; ++i) { // loop through v's neighbors
            uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // u->w with length l
            binfo_t *t = &b->a[w];
            if (w == v0) goto pop_reset;
            if (av[i].del) continue;
            kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
            if (d + l > max_dist) break; // too far
            if (t->s == 0) { // this vertex has never been visited
                kv_push(uint32_t, b->b, w); // save it for revert
                t->p = v, t->s = 1, t->d = d + l;
                t->r = count_out(g, w^1);
                ++n_pending;
            } else { // visited before
                if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
                if (c + 1 > t->c) t->c = c + 1;
                if (d + l < t->d) t->d = d + l; // update dist
            }
            assert(t->r > 0);
            if (--(t->r) == 0) {
                uint32_t x = asg_arc_n(g, w);
                if (x) kv_push(uint32_t, b->S, w);
                else kv_push(uint32_t, b->T, w); // a tip
                --n_pending;
            }
        }
        if (i < nv || b->S.n == 0) goto pop_reset;
    } while (b->S.n > 1 || n_pending);
    asg_bub_backtrack(g, v0, b);
    n_pop = 1 | (uint64_t)b->T.n<<32;
    pop_reset:
    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        binfo_t *t = &b->a[b->b.a[i]];
        t->s = t->c = t->d = 0;
    }
    return n_pop;
}

/// pop bubbles
static int asg_pop_bubble(asg_t *g, int max_dist) {
    uint32_t v, n_vtx = g->n_seq * 2;
    uint64_t n_pop = 0;
    buf_t b;
    if (!g->is_symm) asg_symm(g);
    memset(&b, 0, sizeof(buf_t));
    b.a = hmalloc(binfo_t, n_vtx);
    for (v = 0; v < n_vtx; ++v) {
        uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);
        if (nv < 2 || g->seq[v>>1].del) continue;
        for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
            if (!av[i].del) ++n_arc;
        if (n_arc > 1)
            n_pop += asg_bub_pop1(g, v, max_dist, &b);
    }
    hfree(b.a); hfree(b.S.a); hfree(b.T.a); hfree(b.b.a); hfree(b.e.a);
    if (n_pop) asg_cleanup(g);
    LOG_INFO("[M::%s] popped %d bubbles and trimmed %d tips", __func__, (uint32_t)n_pop, (uint32_t)(n_pop>>32));
    return n_pop;
}

static int asg_cut_internal(asg_t *g, int max_ext) {
    asg64_v a = {0,0,0};
    uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
    for (v = 0; v < n_vtx; ++v) {
        if (g->seq[v>>1].del) continue;
        if (asg_is_utg_end(g, v, 0) != ASG_ET_MULTI_NEI) continue;
        if (asg_extend(g, v, max_ext, &a) != ASG_ET_MULTI_NEI) continue;
        for (i = 0; i < a.n; ++i)
            asg_seq_del(g, (uint32_t)a.a[i]>>1);
        ++cnt;
    }
    hfree(a.a);
    if (cnt > 0) asg_cleanup(g);
    LOG_INFO("[M::%s] cut %d internal sequences", __func__, cnt);
    return cnt;
}

static int asg_cut_biloop(asg_t *g, int max_ext) {
    asg64_v a = {0,0,0};
    uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
    for (v = 0; v < n_vtx; ++v) {
        uint32_t nv, nw, w = UINT32_MAX, x, ov = 0, ox = 0;
        asg_arc_t *av, *aw;
        if (g->seq[v>>1].del) continue;
        if (asg_is_utg_end(g, v, 0) != ASG_ET_MULTI_NEI) continue;
        if (asg_extend(g, v, max_ext, &a) != ASG_ET_MULTI_OUT) continue;
        x = (uint32_t)a.a[a.n - 1] ^ 1;
        nv = asg_arc_n(g, v ^ 1), av = asg_arc_a(g, v ^ 1);
        for (i = 0; i < nv; ++i)
            if (!av[i].del) w = av[i].v ^ 1;
        assert(w != UINT32_MAX);
        nw = asg_arc_n(g, w), aw = asg_arc_a(g, w);
        for (i = 0; i < nw; ++i) { // we are looking for: v->...->x', w->v and w->x
            if (aw[i].del) continue;
            if (aw[i].v == x) ox = aw[i].ol;
            if (aw[i].v == v) ov = aw[i].ol;
        }
        if (ov == 0 && ox == 0) continue;
        if (ov > ox) {
            asg_arc_del(g, w, x, 1);
            asg_arc_del(g, x^1, w^1, 1);
            ++cnt;
        }
    }
    hfree(a.a);
    if (cnt > 0) asg_cleanup(g);
    LOG_INFO("[M::%s] cut %d small bi-loops", __func__, cnt);
    return cnt;
}

/// delete short arcs
static int asg_arc_del_short(asg_t *g, float drop_ratio) {
    uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
    for (v = 0; v < n_vtx; ++v) {
        asg_arc_t *av = asg_arc_a(g, v);
        uint32_t i, thres, nv = asg_arc_n(g, v);
        if (nv < 2) continue;
        thres = (uint32_t)(av[0].ol * drop_ratio + .499);
        for (i = nv - 1; i >= 1 && av[i].ol < thres; --i);
        for (i = i + 1; i < nv; ++i)
            av[i].del = 1, ++n_short;
    }
    if (n_short) {
        asg_cleanup(g);
        asg_symm(g);
    }
    LOG_INFO("[M::%s] removed %d short overlaps", __func__, n_short);
    return n_short;
}

#include "kdq.h"
KDQ_INIT(uint64_t)

#define arc_cnt(g, v) ((uint32_t)(g)->idx[(v)])
#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])

static ma_ug_t *ma_ug_gen(asg_t *g) {
    int32_t *mark;
    uint32_t i, v, n_vtx = g->n_seq * 2;
    kdq_t(uint64_t) *q;
    ma_ug_t *ug;

    ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
    ug->g = hmalloc(asg_t, 1);
    mark = hmalloc(int32_t, n_vtx);

    q = kdq_init(uint64_t);
    for (v = 0; v < n_vtx; ++v) {
        uint32_t w, x, l, start, end, len;
        ma_utg_t *p;
        if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v]) continue;
        mark[v] = 1;
        q->count = 0, start = v, end = v^1, len = 0;
        // forward
        w = v;
        while (1) {
            if (arc_cnt(g, w) != 1) break;
            x = arc_first(g, w).v; // w->x
            if (arc_cnt(g, x^1) != 1) break;
            mark[x] = mark[w^1] = 1;
            l = asg_arc_len(arc_first(g, w));
            kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
            end = x^1, len += l;
            w = x;
            if (x == v) break;
        }
        if (start != (end^1) || kdq_size(q) == 0) { // linear unitig
            l = g->seq[end>>1].len;
            kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
            len += l;
        } else { // circular unitig
            start = end = UINT32_MAX;
            goto add_unitig; // then it is not necessary to do the backward
        }
        // backward
        x = v;
        while (1) { // similar to forward but not the same
            if (arc_cnt(g, x^1) != 1) break;
            w = arc_first(g, x^1).v ^ 1; // w->x
            if (arc_cnt(g, w) != 1) break;
            mark[x] = mark[w^1] = 1;
            l = asg_arc_len(arc_first(g, w));
            kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);
            start = w, len += l;
            x = w;
        }
        add_unitig:
        if (start != UINT32_MAX) mark[start] = mark[end] = 1;
        kv_pushp(ma_utg_t, ug->u, &p);
        p->s = 0, p->start = start, p->end = end, p->len = len, p->n = kdq_size(q), p->circ = (start == UINT32_MAX);
        p->m = p->n;
        kv_roundup32(p->m);
        p->a = (uint64_t*)malloc(8 * p->m);
        for (i = 0; i < kdq_size(q); ++i)
            p->a[i] = kdq_at(q, i);
    }
    kdq_destroy(uint64_t, q);

    // add arcs between unitigs; reusing mark for a different purpose
    for (v = 0; v < n_vtx; ++v) mark[v] = -1;
    for (i = 0; i < ug->u.n; ++i) {
        if (ug->u.a[i].circ) continue;
        mark[ug->u.a[i].start] = i<<1 | 0;
        mark[ug->u.a[i].end] = i<<1 | 1;
    }
    for (i = 0; i < g->n_arc; ++i) {
        asg_arc_t *p = &g->arc[i];
        if (p->del) continue;
        if (mark[p->ul>>32^1] >= 0 && mark[p->v] >= 0) {
            asg_arc_t *q;
            uint32_t u = mark[p->ul>>32^1]^1;
            int l = ug->u.a[u>>1].len - p->ol;
            if (l < 0) l = 1;
            q = asg_arc_pushp(ug->g);
            q->ol = p->ol, q->del = 0;
            q->ul = (uint64_t)u<<32 | l;
            q->v = mark[p->v];
        }
    }
    for (i = 0; i < ug->u.n; ++i)
        asg_seq_set(ug->g, i, ug->u.a[i].len, 0);
    asg_cleanup(ug->g);
    hfree(mark);
    return ug;
}

static void ma_ug_print(const ma_ug_t *ug, const sd_seq_t *seq, int nseq, const ma_sub_t *sub, FILE *fp) {
    uint32_t i, j, l;
    char name[32];
    for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
        ma_utg_t *p = &ug->u.a[i];
        sprintf(name, "utg%.6d%c", i + 1, "lc"[p->circ]);
        fprintf(fp, "S\t%s\t%s\tLN:i:%d\n", name, p->s? p->s : "*", p->len);
        for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
            uint32_t x = p->a[j]>>33;
            if (sub) fprintf(fp, "a\t%s\t%d\tRead-%06d:%d-%d\t%c\t%d\n",
                             name, l, seq[x].id, sub[x].s + 1, sub[x].e, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
            else fprintf(fp, "a\t%s\t%d\tRead-%06d\t%c\t%d\n",
                         name, l, seq[x].id, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
        }
    }
    for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
        uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
        fprintf(fp, "L\tutg%.6d%c\t%c\tutg%.6d%c\t%c\t%dM\tSD:i:%d\n", (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
                (v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
    }
    for (i = 0; i < ug->u.n; ++i) { // summary of unitigs
        uint32_t cnt[2];
        ma_utg_t *u = &ug->u.a[i];
        if (u->start == UINT32_MAX) {
            fprintf(fp, "x\tutg%.6dc\t%d\t%d\n", i + 1, u->len, u->n);
        } else {
            for (j = 0; j < 2; ++j) cnt[j] = asg_arc_n(ug->g, i<<1|j);
            if (sub)
                fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\tRead-%06d:%d-%d\t%c\tRead-%06d:%d-%d\t%c\n", i + 1, u->len, u->n, cnt[1], cnt[0],
                        seq[u->start>>1].id, sub[u->start>>1].s + 1, sub[u->start>>1].e, "+-"[u->start&1],
                        seq[u->end>>1].id, sub[u->end>>1].s + 1, sub[u->end>>1].e, "+-"[u->end&1]);
            else
                fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\tRead-%06d\t%c\tRead-%06d\t%c\n", i + 1, u->len, u->n, cnt[1], cnt[0],
                        seq[u->start>>1].id, "+-"[u->start&1], seq[u->end>>1].id, "+-"[u->end&1]);
        }
    }
}