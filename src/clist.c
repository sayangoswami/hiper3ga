//
// Created by sayan on 6/10/20.
//

#include <hiper3ga.h>
#include <vbyte.h>
#include "clist.h"

static inline block_t* allocate_new(clist_t *self) {
    block_t *b = kdq_pushp(block_t, self);
    b->csize = b->count = 0;
    b->cdata = hmalloc(u_char, CLIST_BLKSZ);
    return b;
}
inline void clist_remove(clist_t *self, int i) {
    if (i < kdq_size(self))
        if (kdq_at(self, i).cdata)
            hfree(kdq_at(self, i).cdata);
}
inline void clist_push_compressed(clist_t *self, const u_char *data, size_t nbytes, size_t count){
    block_t *b = kdq_pushp(block_t, self);
    b->cdata = (void*)data, b->csize = nbytes, b->count = count;
}
inline bool clist_pop_compressed(clist_t *self, u_char **data_p, size_t *nbytes_p, size_t *count_p){
    block_t *b = kdq_shift_block_t(self);
    if (b) {
        *data_p = b->cdata, *nbytes_p = b->csize, *count_p = b->count;
        return true;
    } else return false;
}
inline size_t clist_length(clist_t *self) {
    return kdq_size(self);
}
inline size_t clist_count(clist_t *self) {
    size_t count = 0;
    for (int i = 0; i < kdq_size(self); ++i)
        count += kdq_at(self, i).count;
    return count;
}
inline size_t clist_count_at(clist_t *self, int i) {
    if (i >= kdq_size(self)) return 0;
    else return kdq_at(self, i).count;
}
inline void clist_append(clist_t *one, clist_t *other) {
    assert(one); assert(other);
    if (kdq_size(other)) {
        block_t *b;
        while ((b = kdq_pop_block_t(other)) != NULL)
            kdq_push_block_t(one, *b);
    }
    kdq_destroy_block_t(other);
}

///////////// generic methods /////////////////

inline void clist_push_template(clist_t *self, void *data, size_t count, const size_t itemsize) {
    if (count) {
        block_t *b;
        if (!kdq_size(self)) b = allocate_new(self);
        else b = &kdq_last(self);
        for (int i = 0; i < count; ++i) {
            /// if buffer is full, push to queue and create a new one
            if (GASNETT_PREDICT_FALSE(CLIST_BLKSZ - b->csize < 1.25f * itemsize))
                b = allocate_new(self);
            /// compress and write to buffer.
            b->csize += vbyte_compress_unsorted64((uint64_t*)(data + i * itemsize),
                    b->cdata + b->csize, itemsize / 8);
            b->count++;
        }
    }
}
inline void clist_push_block_template(clist_t *self, void *data, size_t count, void *temp, const size_t itemsize) {
    if (count) {
        size_t csize = vbyte_compress_unsorted64((uint64_t *)data, temp, count * itemsize / 8);
        if (csize > count * itemsize)
            LOG_WARN("Original data size = %zd bytes, compressed size = %zd bytes", count * itemsize, csize);
        block_t *b = kdq_pushp(block_t, self);
        b->count = count, b->csize = csize, b->cdata = hmalloc(u_char, csize);
        memcpy(b->cdata, temp, csize);
    }
}
inline size_t clist_at_template(clist_t *self, int i, void *buffer, size_t bufsz, const size_t itemsize) {
    if (i >= kdq_size(self)) return 0;
    else {
        assert(buffer && bufsz);
        block_t *b = &kdq_at(self, i);
        if (bufsz < b->count)
            LOG_ERR("Buffer is not large enough.");
        size_t csize = vbyte_uncompress_unsorted64(b->cdata, (uint64_t*)(buffer), b->count * itemsize / 8);
        assert(csize == b->csize);
        return b->count;
    }
}
inline size_t clist_pop_template(clist_t *self, void *buffer, size_t bufsz, const size_t itemsize) {
    block_t *b = kdq_shift_block_t(self);
    if (b) {
        if (bufsz < b->count)
            LOG_ERR("Buffer is not large enough.");
        size_t csize = vbyte_uncompress_unsorted64(b->cdata, (uint64_t*)(buffer), b->count * itemsize/8);
        assert(csize == b->csize);
        hfree(b->cdata);
        return b->count;
    } else return 0;
}
inline double clist_compfactor_template(clist_t *self, const size_t itemsize){
    size_t count = 0, csize = 0;
    for (int i = 0; i < kdq_size(self); ++i)
        count += kdq_at(self, i).count, csize += kdq_at(self, i).csize;
    return csize * 1.0/(1.0 * count * itemsize);
}