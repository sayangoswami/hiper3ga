//
// Created by sayan on 9/13/19.
//

#ifndef HIPER3GA_HBITVEC_H
#define HIPER3GA_HBITVEC_H

hbitvec_t *hbitvec_new(size_t size);

bool hbitvec_is(void *self);

void hbitvec_destroy(hbitvec_t **self_p);

bool hbitvec_set(hbitvec_t *self, size_t i);

bool hbitvec_get(hbitvec_t *self, size_t i);

uint32_t* hbitvec_raw(hbitvec_t *self, size_t *size_p);

void hbitvec_clear_all(hbitvec_t *self);

#endif //HIPER3GA_HBITVEC_H
