#ifndef ABPOA_SIMD_H
#define ABPOA_SIMD_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const int reg_n, bits_n, log_num, num_of_value, size;
    int inf_min; // based on penalty of mismatch and GAP_OE1
} SIMD_para_t;

#define SIMDShiftOneNi8  1
#define SIMDShiftOneNi16 2
#define SIMDShiftOneNi32 4
#define SIMDShiftOneNi64 8

#if __AVX512BW__
#define SIMDTotalBytes 64
#elif defined(__AVX2__)
#define SIMDTotalBytes 32
#else
#define SIMDTotalBytes 16
#endif


abpoa_simd_matrix_t *abpoa_init_simd_matrix(void);
void abpoa_free_simd_matrix(abpoa_simd_matrix_t *abm);
void simd_output_pre_nodes(int *pre_index, int pre_n, int dp_i, int dp_j, int cur_op, int verbose);
int simd_abpoa_realloc(abpoa_t *ab, int gn, int qlen, abpoa_para_t *abpt, SIMD_para_t sp);

#ifdef __cplusplus
}
#endif

#endif
