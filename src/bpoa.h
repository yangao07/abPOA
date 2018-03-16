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
    int align_mode; // 0: global, 1: local, 2: extend
    // available SIMD instruction
    int simd_flag;
} bpoa_para_t;

bpoa_para_t *bpoa_para_init(void);
void bpoa_para_free(bpoa_para_t *bpt);

// int bpoa_main(const char *seq_fn, bpoa_para_t *bpt) { TODO
int bpoa_main(int seq_n, char (*seq)[100], bpoa_para_t *bpt);

#ifdef __cplusplus
}
#endif

#endif
