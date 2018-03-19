#ifndef POA_H
#define POA_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    // score matrix
    int m; int *mat;
    int match, mismatch, gap_open, gap_ext; int inf_min;
    int bw; // band width
    int zdrop, end_bonus; // from minimap2
    // alignment mode
    int8_t align_mode:2, use_ada:2; // mode: 0: global, 1: local, 2: extend
    // available SIMD instruction
    int simd_flag;
} abpoa_para_t;

abpoa_para_t *abpoa_para_init(void);
void abpoa_para_free(abpoa_para_t *abpt);

// int abpoa_main(const char *seq_fn, abpoa_para_t *abpt) { TODO
int abpoa_main(int seq_n, char (*seq)[100], abpoa_para_t *abpt);

#ifdef __cplusplus
}
#endif

#endif
