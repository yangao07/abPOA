#ifndef ABPOA_SIMD_INTERNAL_H
#define ABPOA_SIMD_INTERNAL_H

#include "simd_instruction.h"

struct abpoa_simd_matrix_t {
    SIMDi *s_mem; uint64_t s_msize;
    int *dp_beg, *dp_end, *dp_beg_sn, *dp_end_sn, rang_m;
};

#endif
