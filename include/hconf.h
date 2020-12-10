//
// Created by sayan on 9/13/19.
//

#ifndef HIPER3GA_HCONF_H
#define HIPER3GA_HCONF_H

extern bool bi_dir;
extern int min_span;
extern int min_match;
extern int min_dp;
extern float min_iden;

extern int max_hang;
extern int min_ovlp;
extern float int_frac;

extern int gap_fuzz;
extern int n_rounds;
extern int bub_dist;
extern int max_ext;
extern float min_ovlp_drop_ratio;
extern float max_ovlp_drop_ratio;
extern float final_ovlp_drop_ratio;

#define ENCODER_VB 0x1
#define ENCODER_VBX 0x2
#define ENCODER_VBZ 0x3
#define ENCODER_VS 0x4
#define ENCODER_P4 0x5
#define ENCODER_P4Z 0x6
#define ENCODER_BP 0x7
#define ENCODER_BPZ 0x8
//// delta
#define ENCODER_VBD 0x10
#define ENCODER_P4D 0x11
#define ENCODER_BPD 0x12
#define ENCODER_BPF 0x13
#define ENCODER_EF 0x14

#define ASM_NPTN 256

hconf_t *hconf_create(int argc, char **argv, hctx_t *ctx);

void hconf_destroy(hconf_t **self_p);

bool hconf_is(void *self);

char* hconf_input(hconf_t *self);

char* hconf_output(hconf_t *self);

size_t hconf_max_bases_per_batch(hconf_t *self);

char* hconf_tempdir(hconf_t *self);

int hconf_nthreads(hconf_t *self);

int hconf_niothreads(hconf_t *self);

int hconf_window_len(hconf_t *self);

int hconf_kmer_len(hconf_t *self);

bool hconf_sync_requests();

#endif //HIPER3GA_HCONF_H
