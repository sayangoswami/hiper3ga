//
// Created by sayan on 6/11/20.
//

#ifndef HIPER3GA_PARTITIONER_H
#define HIPER3GA_PARTITIONER_H

#include "hiper3ga.h"
#include "clist.h"
#include "ksort.h"
#include "kthread.h"

#define HPTN_INIT(name, type, key, keysize, NBUCKETS, MAX_TBUF_SZ)                                      \
    typedef struct {                                                                                    \
        type *buffer, *temp;                                                                            \
        size_t bufcount, tempcount;                                                                     \
        clist_t **listout;                                                                              \
    } _##name##_tbuf_t;                                                                                 \
    typedef struct {                                                                                    \
        runconf_t *runconf;                                                                             \
        clist_t *listin;                                                                                \
        _##name##_tbuf_t *tbuf;                                                                         \
    } _##name##_shared_data_t;                                                                          \
    static inline __attribute__((always_inline)) int _##name##_ptnkey(type a)                           \
                                                                    { return key(a) % NBUCKETS; }       \
    KRADIX_SORT_INIT(name, type, _##name##_ptnkey, keysize)                                             \
    static void _##name##_partitioner_fn(void *s, long ignore, int tid) {                               \
        _##name##_shared_data_t *sd = (_##name##_shared_data_t*)s;                                      \
        _##name##_tbuf_t *tbuf = &sd->tbuf[tid];                                                        \
        tbuf->buffer = hmalloc(type, MAX_TBUF_SZ), tbuf->bufcount = 0;                                  \
        tbuf->temp = NULL, tbuf->tempcount = 0;                                                         \
        tbuf->listout = hmalloc(clist_t*, NBUCKETS);                                                    \
        for (int i = 0; i < NBUCKETS; ++i) tbuf->listout[i] = clist_create();                           \
        uint ptnsz[NBUCKETS];                                                                           \
        for (int t = tid; t < clist_length(sd->listin); t += sd->runconf->NTHREADS) {                   \
            size_t c = clist_count_at(sd->listin, t);                                                   \
            if (tbuf->bufcount + c > MAX_TBUF_SZ) {                                                     \
                /** partition and save into individual lists */\
                radix_sort_##name(tbuf->buffer, tbuf->buffer + tbuf->bufcount);                         \
                memset(ptnsz, 0, 4*NBUCKETS);                                                           \
                for (uint j = 0; j < tbuf->bufcount; ++j)                                               \
                    ptnsz[_##name##_ptnkey(tbuf->buffer[j])]++;                                                   \
                for (uint j = 0, o = 0; j < NBUCKETS; o += ptnsz[j], ++j) {                             \
                    if (tbuf->tempcount < 2*ptnsz[j]) {                                                   \
                        tbuf->temp = hrealloc(tbuf->temp, type, 2*ptnsz[j]);                            \
                        tbuf->tempcount = 2*ptnsz[j];                                                   \
                    }                                                                                   \
                    /** compress and push compressed */\
                    clist_push_block(tbuf->listout[j], type, tbuf->buffer + o, ptnsz[j], tbuf->temp);   \
                }                                                                                       \
                tbuf->bufcount = 0;                                                                     \
            }                                                                                           \
            clist_at(sd->listin, type, t, tbuf->buffer + tbuf->bufcount, c);                            \
            tbuf->bufcount += c;                                                                        \
            clist_remove(sd->listin, t);                                                                \
        }                                                                                               \
        if (tbuf->bufcount) {                                                                           \
            /** partition and save into individual lists */\
            radix_sort_##name(tbuf->buffer, tbuf->buffer + tbuf->bufcount);                             \
            memset(ptnsz, 0, 4*NBUCKETS);                                                               \
            for (uint j = 0; j < tbuf->bufcount; ++j)                                                   \
                ptnsz[_##name##_ptnkey(tbuf->buffer[j])]++;                                                       \
            for (uint j = 0, o = 0; j < NBUCKETS; o += ptnsz[j], ++j) {                                 \
                if (tbuf->tempcount < 2*ptnsz[j]) {                                                       \
                    tbuf->temp = hrealloc(tbuf->temp, type, 2*ptnsz[j]);                                \
                    tbuf->tempcount = 2*ptnsz[j];                                                       \
                }                                                                                       \
                /** compress and push compressed */\
                clist_push_block(tbuf->listout[j], type, tbuf->buffer + o, ptnsz[j], tbuf->temp);       \
            }                                                                                           \
        }                                                                                               \
        hfree(tbuf->buffer); hfree(tbuf->temp);                                                         \
    }                                                                                                   \
    static clist_t ** _##name##_partition(clist_t *listin, runconf_t *runconf) {                        \
        LOG_INFO("Partitioning..");                                                                     \
        _##name##_tbuf_t tbuf[runconf->NTHREADS];                                                       \
        _##name##_shared_data_t sd = {                                                                  \
                .runconf = runconf, .tbuf = tbuf, .listin = listin,                                     \
        };                                                                                              \
        kt_forpool(runconf->tpool, _##name##_partitioner_fn, &sd, runconf->NTHREADS);                   \
        clist_t ** listout = hmalloc(clist_t*, NBUCKETS);                                               \
        for (int i = 0; i < NBUCKETS; ++i) listout[i] = clist_create();                                 \
        u_char *data; size_t size, count;                                                               \
        for (int i = 0; i < runconf->NTHREADS; ++i) {                                                   \
            for (int j = 0; j < NBUCKETS; ++j) {                                                        \
                while (clist_pop_compressed(tbuf[i].listout[j], &data, &size, &count))                  \
                    clist_push_compressed(listout[j], data, size, count);                               \
                clist_destroy(tbuf[i].listout[j]);                                                      \
            }                                                                                           \
            hfree(tbuf[i].listout);                                                                     \
        }                                                                                               \
        return listout;                                                                                 \
    }

#define partition(name, listin, runconf) _##name##_partition(listin, runconf)

#endif //HIPER3GA_PARTITIONER_H
