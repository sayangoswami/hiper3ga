//
// Created by sayan on 6/10/20.
//

#ifndef HIPER3GA_CLIST_H
#define HIPER3GA_CLIST_H

#include "kdq.h"
#include "hfile.h"

#define CLIST_BLKSZ (4u MiB)

typedef struct {
    void *cdata;
    uint32_t count, csize;
} block_t;
KDQ_INIT(block_t)
typedef kdq_t(block_t) clist_t;

void clist_remove(clist_t *self, int i);
void clist_push_compressed(clist_t *self, const u_char *data, size_t nbytes, size_t count);
bool clist_pop_compressed(clist_t *self, u_char **data_p, size_t *nbytes_p, size_t *count_p);
size_t clist_length(clist_t *self);
size_t clist_count(clist_t *self);
size_t clist_count_at(clist_t *self, int i);
void clist_append(clist_t *one, clist_t *other);

void clist_push_template(clist_t *self, void *data, size_t count, size_t itemsize);
void clist_push_block_template(clist_t *self, void *data, size_t count, void *temp, size_t itemsize);
size_t clist_at_template(clist_t *self, int i, void *buffer, size_t bufsz, size_t itemsize);
size_t clist_pop_template(clist_t *self, void *buffer, size_t bufsz, size_t itemsize);
double clist_compfactor_template(clist_t *self, size_t itemsize);

#define clist_create()      \
        kdq_init_block_t()

#define clist_destroy(p_list)   \
        kdq_destroy_block_t(p_list)

#define clist_push(p_list, type, data, count)       \
        clist_push_template(p_list, data, count, sizeof(type))

#define clist_push_block(p_list, type, data, count, temp)       \
        clist_push_block_template(p_list, data, count, temp, sizeof(type))

#define clist_at(p_list, type, i, buffer, bufsz)  \
        clist_at_template(p_list, i, buffer, bufsz, sizeof(type))

#define clist_pop(p_list, type, buffer, bufsz) \
        clist_pop_template(p_list, buffer, bufsz, sizeof(type))

#define clist_compfactor(p_list, type)              \
        clist_compfactor_template(p_list, sizeof(type))


#endif //HIPER3GA_CLIST_H
