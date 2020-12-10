//
// Created by sayan on 6/5/20.
//

#include "partitioner.h"
#include "hidx.h"

static int NNODES;

#define u128L1_key(a) ((a).x)
HPTN_INIT(u128L1, u128_t, u128L1_key, 4, NNODES, 16 MiB)
clist_t ** u128L1_partition(clist_t *listin, runconf_t *runconf) {
    NNODES = runconf->NNODES;
    return partition(u128L1, listin, runconf);
}

#define u128L2_key(a) ((a).x /  NNODES)
HPTN_INIT(u128L2, u128_t, u128L2_key, 4, NPTN, 16 MiB)
clist_t ** u128L2_partition(clist_t *listin, runconf_t *runconf) {
    NNODES = runconf->NNODES;
    return partition(u128L2, listin, runconf);
}

#define hitL1_key(a) ((uint32_t)((a).qns>>32u))
HPTN_INIT(hitL1, ma_hit_t, hitL1_key, 4, NNODES, 16 MiB)

clist_t ** hitL1_partition(clist_t *listin, runconf_t *runconf) {
    NNODES = runconf->NNODES;
    return partition(hitL1, listin, runconf);
}

#define hitL2_key(a) ((uint32_t)((a).qns>>32u) / NNODES)
HPTN_INIT(hitL2, ma_hit_t, hitL2_key, 4, ASM_NPTN, 8 MiB)

clist_t ** hitL2_partition(clist_t *listin, runconf_t *runconf) {
    NNODES = runconf->NNODES;
    return partition(hitL2, listin, runconf);
}