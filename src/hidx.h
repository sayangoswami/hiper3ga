//
// Created by sayan on 5/31/20.
//

#ifndef HIPER3GA_HIDX_H
#define HIPER3GA_HIDX_H

#include "hiper3ga.h"
#include "clist.h"
#include "metadata.h"

#define PBITS 14u
#define NPTN (1u<<PBITS)
#define PTNMASK (NPTN-1)

void * idx_create(clist_t** lists, readinfo_t *readinfo, hctx_t *ctx, runconf_t *runconf);
void idx_destroy(void **self_p);
void idx_setup_request(void *p_idx, int client_id, u128_t *mm, size_t nmm);
void idx_send_requests(void *p_idx);
void idx_setup_response(void *p_idx, int client_id);
unsigned idx_nreads(void *p_idx, int client_id);
u128_t * idx_get_mm(void *p_idx, int client_id, int i, uint *p_nmm, uint *p_rid);
const uint64_t * idx_get_vals(void *p_idx, int client_id, uint64_t key, uint16_t *n_vals);
void idx_wait(void *p_idx);

#define mm_rid(a) (((a).y)>>32u)
#define mm_pos(a) ((uint32_t)((a).y))

#endif //HIPER3GA_HIDX_H
