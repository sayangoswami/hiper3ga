//
// Created by sayan on 9/13/19.
//

#include "hclasses.h"

/** Sets bit at position V of the value at address P and returns true if it was previously set */
#define atomic_set_bit(P, V)  (__atomic_fetch_or((P), 0x80000000 >> (V), __ATOMIC_RELAXED) & (0x80000000 >> (V)))

/** Returns the bit at position V of the value at address P */
#define atomic_get_bit(P, V) (__atomic_load_n((P), __ATOMIC_RELAXED) & (0x80000000 >> (V)))

/** Clears the bit at position V of the value at address P */
#define atomic_clear_bit(P, V) (__atomic_fetch_and((P), ~(0x80000000 >> (V)), __ATOMIC_RELAXED))

/** Compile read-write barrier */
#define barrier() asm volatile("": : :"memory")

/** Pause instruction to prevent excess processor bus usage */
#define cpu_relax() asm volatile("pause\n": : :"memory")

static const char *TAGSTR = "BTVC";

struct _hbitvec_t {
    uint32_t tag;
    uint32_t *array;
    size_t size;
    size_t array_size;
};

hbitvec_t *hbitvec_new(size_t size) {
    hbitvec_t *self = hmalloc(hbitvec_t, 1);
    self->tag = *((uint32_t *)TAGSTR);
    self->size = size;
    self->array_size = size/32 + (size%32 != 0);
    self->array = hmalloc(uint32_t, self->array_size);
    return self;
}

bool hbitvec_is(void *self) {
    assert(self);
    return ((hbitvec_t*)self)->tag == *((uint32_t *)TAGSTR);
}

void hbitvec_destroy(hbitvec_t **self_p) {
    if (self_p && *self_p) {
        hbitvec_t *self = *self_p;
        assert(hbitvec_is(self));
        hfree(self->array);
        hfree(self);
        *self_p = NULL;
    }
}

bool hbitvec_set(hbitvec_t *self, size_t i) {
    assert(hbitvec_is(self));
    assert(i < self->size);
    size_t offset = i / 32;
    size_t bitpos = i % 32;
    return atomic_set_bit(self->array + offset, bitpos);
}

bool hbitvec_get(hbitvec_t *self, size_t i) {
    assert(hbitvec_is(self));
    assert(i < self->size);
    size_t offset = i / 32;
    size_t bitpos = i % 32;
    return atomic_get_bit(self->array + offset, bitpos);
}

uint32_t* hbitvec_raw(hbitvec_t *self, size_t *size_p) {
    assert(hbitvec_is(self));
    assert(size_p);
    *size_p = self->size;
    return self->array;
}

void hbitvec_clear_all(hbitvec_t *self) {
    assert(hbitvec_is(self));
    memset(self->array, 0, self->array_size * sizeof(uint32_t));
}
