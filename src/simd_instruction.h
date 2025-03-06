// A header file to get you set going with Intel SIMD instrinsic programming. 
// <immintrin.h> is inlucded for SSE2, SSE41, AVX2 and AVX512F, AVX512BW
// SSE4.1: floor and blend is available)
// AVX2: double speed

// do not support AVX512F/AVX512BW 12/20/2021 - Yan Gao
// AVX512F: quardruple speed
// AVX512BW: byte and word operation

#include <stdlib.h>
#include <errno.h>

#pragma once
#ifndef SIMD_INSTRUCTION_H
#define SIMD_INSTRUCTION_H


#ifndef USE_SIMDE
#include <immintrin.h>
#else // use SIMDE
#ifdef __AVX512BW__
#include "simde/simde/x86/avx512.h"
#else
#ifdef __AVX2__
#include "simde/simde/x86/avx2.h"
#else
#ifdef __SSE4_1__
#include "simde/simde/x86/sse4.1.h"
#else
#include "simde/simde/x86/sse2.h"
#endif // end of sse41
#endif // end of AVX2
#endif // end of 512F
#endif // end of USE_SIMDE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define SIMD_SSE      0x1
#define SIMD_SSE2     0x2
#define SIMD_SSE3     0x4
#define SIMD_SSSE3    0x8
#define SIMD_SSE41    0x10
#define SIMD_SSE42    0x20
#define SIMD_AVX      0x40
#define SIMD_AVX2     0x80
#define SIMD_AVX512F  0x100
#define SIMD_AVX512BW 0x200

// #define SIMDFree(x) _mm_free(x)
// posix_memalign and free
#define SIMDFree(x) free(x)

// Shift, Blend, ... for 8/16 and 32/64
#ifdef __AVX512BW__ // start of AVX512BW

typedef __m512 SIMDf;
typedef __m512i SIMDi;

#define SIMDStore(x,y) _mm512_store_ps(x,y)
#define SIMDStorei(x,y) _mm512_store_si512(x,y)
#define SIMDLoad(x) _mm512_load_ps(x)
#define SIMDLoadi(x) _mm512_load_si512(x)
#define SIMDZero _mm512_setzero_si512()
#define SIMDSetZero() _mm512_setzero_ps()
#define SIMDSetZeroi() _mm512_setzero_si512()
#define SIMDSetOne(x) _mm512_set1_ps(x)
#define SIMDSetOnei8(x) _mm512_set1_epi8(x)
#define SIMDSetOnei16(x) _mm512_set1_epi16(x)
#define SIMDSetOnei32(x) _mm512_set1_epi32(x)
#define SIMDSetOnei64(x) _mm512_set1_epi64(x)
#define SIMDAdd(x,y) _mm512_add_ps(x,y)
#define SIMDAddi8(x,y) _mm512_add_epi8(x,y)
#define SIMDAddi16(x,y) _mm512_add_epi16(x,y)
#define SIMDAddi32(x,y) _mm512_add_epi32(x,y)
#define SIMDAddi64(x,y) _mm512_add_epi64(x,y)
#define SIMDSub(x,y) _mm512_sub_ps(x,y)
#define SIMDSubi8(x,y) _mm512_sub_epi8(x,y)
#define SIMDSubi16(x,y) _mm512_sub_epi16(x,y)
#define SIMDSubi32(x,y) _mm512_sub_epi32(x,y)
#define SIMDSubi64(x,y) _mm512_sub_epi64(x,y)
#define SIMDMul(x,y) _mm512_mul_ps(x,y)
#define SIMDMuli32(x,y) _mm512_mul_epi32(x,y)
#define SIMDAnd(x,y) _mm512_and_ps(x,y)
#define SIMDAndi(x,y) _mm512_and_si512(x,y)
#define SIMDAndNot(x,y) _mm512_andnot_ps(x,y)
#define SIMDAndNoti(x,y) _mm512_andnot_si512(x,y)
#define SIMDOr(x,y) _mm512_or_ps(x,y)
#define SIMDOri(x,y) _mm512_or_si512(x,y)
#define SIMDShiftLeft(x,n) \
    (n) < 16 ? \
    _mm512_alignr_epi8(x, _mm512_shuffle_i64x2(_mm512_shuffle_i64x2(x, SIMDZero, _MM_SHUFFLE(0,0,1,0)), x, _MM_SHUFFLE(2,1,0,2)), (16-(n))) : \
    ((n) < 32 ? \
    _mm512_alignr_epi8(_mm512_shuffle_i64x2(_mm512_shuffle_i64x2(x, SIMDZero, _MM_SHUFFLE(0,0,1,0)), x, _MM_SHUFFLE(2,1,0,2)), _mm512_shuffle_i64x2(SIMDZero, x, _MM_SHUFFLE(1,0,0,0)), (32-(n))) : \
    ((n) < 48 ? \
    _mm512_alignr_epi8(_mm512_shuffle_i64x2(SIMDZero, x, _MM_SHUFFLE(1,0,0,0)), _mm512_shuffle_i64x2(SIMDZero, _mm512_shuffle_i64x2(SIMDZero, x, _MM_SHUFFLE(1,0,0,0)), _MM_SHUFFLE(2,0,0,0)), (48-(n))) : \
    _mm512_bslli_epi128(_mm512_shuffle_i64x2(SIMDZero,  _mm512_shuffle_i64x2(SIMDZero, x, _MM_SHUFFLE(1,0,0,0)), _MM_SHUFFLE(2,0,0,0)), ((n)-48))))

#define SIMDShiftRight(x,n) \
    (n) < 16 ? \
    _mm512_alignr_epi8(_mm512_shuffle_i64x2( _mm512_shuffle_i64x2(SIMDZero, x, _MM_SHUFFLE(3,2,0,0)), x, _MM_SHUFFLE(0,3,2,1)), x, (n)) : \
    ((n) < 32 ? \
    _mm512_alignr_epi8(_mm512_shuffle_i64x2(x, SIMDZero, _MM_SHUFFLE(0,0,3,2)), _mm512_shuffle_i64x2(_mm512_shuffle_i64x2(SIMDZero, x, _MM_SHUFFLE(3,2,0,0)), x, _MM_SHUFFLE(0,3,2,1)), ((n)-16)) : \
    ((n) < 48 ? \
    _mm512_alignr_epi8(_mm512_shuffle_i64x2(_mm512_shuffle_i64x2(x, SIMDZero, _MM_SHUFFLE(0,0,3,2)), SIMDZero, _MM_SHUFFLE(0,0,2,1)), _mm512_shuffle_i64x2(x, SIMDZero, _MM_SHUFFLE(0,0,3,2)), ((n)-32)) : \
    _mm512_bsrli_epi128(_mm512_shuffle_i64x2(_mm512_shuffle_i64x2(x, SIMDZero, _MM_SHUFFLE(0,0,3,2)), SIMDZero, _MM_SHUFFLE(0,0,2,1)), ((n)-48))))

#define SIMDShiftLeftOnei16(x,y) _mm512_slli_epi16(x,y)
#define SIMDShiftLeftOnei32(x,y) _mm512_slli_epi32(x,y)
#define SIMDShiftLeftOnei64(x,y) _mm512_slli_epi64(x,y)
#define SIMDShiftRightOnei16(x,y) _mm512_srli_epi16(x,y)
#define SIMDShiftRightOnei32(x,y) _mm512_srli_epi32(x,y)
#define SIMDShiftRightOnei64(x,y) _mm512_srli_epi64(x,y)
#define SIMDEqualM(x,y) _mm512_cmpeq_ps_mask(x,y)
#define SIMDEquali8M(x,y) _mm512_cmpeq_epi8_mask(x,y)
#define SIMDEquali16M(x,y) _mm512_cmpeq_epi16_mask(x,y)
#define SIMDEquali32M(x,y) _mm512_cmpeq_epi32_mask(x,y)
#define SIMDEquali64M(x,y) _mm512_cmpeq_epi64_mask(x,y)
#define SIMDNotEqualM(x,y) _mm512_cmpneq_ps_mask(x,y)
#define SIMDNotEquali8M(x,y) _mm512_cmpneq_epi8_mask(x,y)
#define SIMDNotEquali16M(x,y) _mm512_cmpneq_epi16_mask(x,y)
#define SIMDNotEquali32M(x,y) _mm512_cmpneq_epi32_mask(x,y)
#define SIMDNotEquali64M(x,y) _mm512_cmpneq_epi64_mask(x,y)
#define SIMDGreaterThani8M(x,y) _mm512_cmpgt_epi8_mask(x,y)
#define SIMDGreaterThani16M(x,y) _mm512_cmpgt_epi16_mask(x,y)
#define SIMDGreaterThani32M(x,y) _mm512_cmpgt_epi32_mask(x,y)
#define SIMDGreaterThani64M(x,y) _mm512_cmpgt_epi64_mask(x,y)
#define SIMDGreaterThanOrEquali8M(x,y) _mm512_cmpge_epi8_mask(x,y)
#define SIMDGreaterThanOrEquali16M(x,y) _mm512_cmpge_epi16_mask(x,y)
#define SIMDGreaterThanOrEquali32M(x,y) _mm512_cmpge_epi32_mask(x,y)
#define SIMDGreaterThanOrEquali64M(x,y) _mm512_cmpge_epi64_mask(x,y)
#define SIMDLessThanM(x,y) _mm512_cmplt_ps_mask(x,y)
#define SIMDLessThani8M(x,y) _mm512_cmplt_epi8_mask(x,y)
#define SIMDLessThani16M(x,y) _mm512_cmplt_epi16_mask(x,y)
#define SIMDLessThani32M(x,y) _mm512_cmplt_epi32_mask(x,y)
#define SIMDLessThani64M(x,y) _mm512_cmplt_epi64_mask(x,y)
#define SIMDLessThanOrEqualM(x,y) _mm512_cmple_ps_mask(x,y)
#define SIMDLessThanOrEquali8M(x,y) _mm512_cmple_epi8_mask(x,y)
#define SIMDLessThanOrEquali16M(x,y) _mm512_cmple_epi16_mask(x,y)
#define SIMDLessThanOrEquali32M(x,y) _mm512_cmple_epi32_mask(x,y)
#define SIMDLessThanOrEquali64M(x,y) _mm512_cmple_epi64_mask(x,y)
#define SIMDMax(x,y) _mm512_max_ps(x,y)
#define SIMDMaxi8(x,y) _mm512_max_epi8(x,y)
#define SIMDMaxi16(x,y) _mm512_max_epi16(x,y)
#define SIMDMaxi32(x,y) _mm512_max_epi32(x,y)
#define SIMDMaxi64(x,y) _mm512_max_epi64(x,y)
#define SIMDMin(x,y) _mm512_min_ps(x,y)
#define SIMDMini8(x,y) _mm512_min_epi8(x,y)
#define SIMDMini16(x,y) _mm512_min_epi16(x,y)
#define SIMDMini32(x,y) _mm512_min_epi32(x,y)
#define SIMDMini64(x,y) _mm512_min_epi64(x,y)

#define SIMDBlend(x,y,z) _mm512_mask_blend_ps(z, x, y)
#define SIMDBlendi8(x,y,z) _mm512_mask_blend_epi8(z, x, y)
#define SIMDBlendi16(x,y,z) _mm512_mask_blend_epi16(z, x, y)
#define SIMDBlendi32(x,y,z) _mm512_mask_blend_epi32(z, x, y)
#define SIMDBlendi64(x,y,z) _mm512_mask_blend_epi64(z, x, y)

// with AVX512BW
#define Maski8 __mmask64
#define Maski16 __mmask32
#define Maski32 __mmask16
#define Maski64 __mmask8
/* x = a == b ? c : d */ 
#define SIMDSetIfEquali8(x,a,b,c,d)  { x = SIMDBlendi8(d, c, SIMDEquali8M(a,b)); } 
#define SIMDSetIfEquali16(x,a,b,c,d) { x = SIMDBlendi16(d, c, SIMDEquali16M(a,b)); } 
#define SIMDSetIfEquali32(x,a,b,c,d) { x = SIMDBlendi32(d, c, SIMDEquali32M(a,b)); } 
#define SIMDSetIfEquali64(x,a,b,c,d) { x = SIMDBlendi64(d, c, SIMDEquali64M(a,b)); } 
/* x = a > b ? c : d */
#define SIMDSetIfGreateri8(x,a,b,c,d)  { x = SIMDBlendi8(d, c, SIMDGreaterThani8M(a,b)); } 
#define SIMDSetIfGreateri16(x,a,b,c,d) { x = SIMDBlendi16(d, c, SIMDGreaterThani16M(a,b)); } 
#define SIMDSetIfGreateri32(x,a,b,c,d) { x = SIMDBlendi32(d, c, SIMDGreaterThani32M(a,b)); } 
#define SIMDSetIfGreateri64(x,a,b,c,d) { x = SIMDBlendi32(d, c, SIMDGreaterThani64M(a,b)); } 
/* x = a < b ? c : d */
#define SIMDSetIfLessi8(x,a,b,c,d)  { x = SIMDBlendVi8(d, c, SIMDGreaterThani8M(b,a)); } 
#define SIMDSetIfLessi16(x,a,b,c,d) { x = SIMDBlendi16(d, c, SIMDGreaterThani16M(b,a)); } 
#define SIMDSetIfLessi32(x,a,b,c,d) { x = SIMDBlendi32(d, c, SIMDGreaterThani32M(b,a)); } 
#define SIMDSetIfLessi64(x,a,b,c,d) { x = SIMDBlendi64(d, c, SIMDGreaterThani64M(b,a)); } 

/* x = a > b ? c : d, y = a > b ? a : b */
#define SIMDGetIfGreateri8(x,y,a,b,c,d)  { Maski8  cmp = SIMDGreaterThani8M(a,b);  x = SIMDBlendi8(d, c, cmp);  y = SIMDBlendi8(b, a, cmp); } 
#define SIMDGetIfGreateri16(x,y,a,b,c,d) { Maski16 cmp = SIMDGreaterThani16M(a,b); x = SIMDBlendi16(d, c, cmp); y = SIMDBlendi16(b, a, cmp); } 
#define SIMDGetIfGreateri32(x,y,a,b,c,d) { Maski32 cmp = SIMDGreaterThani32M(a,b); x = SIMDBlendi32(d, c, cmp); y = SIMDBlendi32(b, a, cmp); } 
#define SIMDGetIfGreateri64(x,y,a,b,c,d) { Maski64 cmp = SIMDGreaterThani64M(a,b); x = SIMDBlendi64(d, c, cmp); y = SIMDBlendi64(b, a, cmp); } 
/* x = a < b ? c : d, y = a < b ? a : b */
#define SIMDGetIfLessi8(x,y,a,b,c,d)  { Maski8  cmp = SIMDGreaterThani8M(b,a);  x = SIMDBlendi8(d, c, cmp);  y = SIMDBlendi8(b, a, cmp); } 
#define SIMDGetIfLessi16(x,y,a,b,c,d) { Maski16 cmp = SIMDGreaterThani16M(b,a); x = SIMDBlendi16(d, c, cmp); y = SIMDBlendi16(b, a, cmp); } 
#define SIMDGetIfLessi32(x,y,a,b,c,d) { Maski32 cmp = SIMDGreaterThani32M(b,a); x = SIMDBlendi32(d, c, cmp); y = SIMDBlendi32(b, a, cmp); } 
#define SIMDGetIfLessi64(x,y,a,b,c,d) { Maski64 cmp = SIMDGreaterThani64M(b,a); x = SIMDBlendi64(d, c, cmp); y = SIMDBlendi64(b, a, cmp); } 

// AVX512F has no  following instructions (AVX512BW HAS), so AVX512F (without AVX512BW) is not working for 8/16 bits tasks
// addi8/16, subi8/16, alignri8, bslli_epi128, bslrli_epi128,
// comeqi8/16, cmpneqi8/16, cmpgti8/16, cmpgei8/16, cmplti8/16, cmplei8
// maxi8/16, blendi8/i16, slli_epi16,srli_epi16 
// end of AVX512BW
#else  // AVX2 SSE4.1 SSE2
#ifdef __AVX2__

// start of AVX2
// m256 will be our base type
typedef __m256 SIMDf;  //for floats
typedef __m256i SIMDi; //for integers

//intrinsic functions
#define SIMDStore(x,y) _mm256_store_ps(x,y)
#define SIMDLoad(x) _mm256_load_ps(x)
#define SIMDStorei(x,y) _mm256_store_si256(x,y)
#define SIMDLoadi(x) _mm256_load_si256(x)
#define SIMDSet(x,y,z,w,a,b,c,d) _mm256_set_ps(x,y,z,w,a,b,c,d)
#define SIMDSeti8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32) __mm256_set_epi8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32)
#define SIMDSeti16(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16) __mm256_set_epi16(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16)
#define SIMDSeti32(x1,x2,x3,x4,x5,x6,x7,x8) __mm256_set_epi32(x1,x2,x3,x4,x5,x6,x7,x8)
#define SIMDSeti64(x1,x2,x3,x4) __mm256_set_epi64x(x1,x2,x3,x4)
#define SIMDSeti128(x,y) __mm256_set_m128(x,y)
#define SIMDSetZero() _mm256_setzero_ps()
#define SIMDSetZeroi() _mm256_setzero_si256()
#define SIMDSetOne(x) _mm256_set1_ps(x)
#define SIMDSetOnei8(x) _mm256_set1_epi8(x)
#define SIMDSetOnei16(x) _mm256_set1_epi16(x)
#define SIMDSetOnei32(x) _mm256_set1_epi32(x)
#define SIMDSetOnei64(x) _mm256_set1_epi64x(x)
#define SIMDAdd(x,y) _mm256_add_ps(x,y)
#define SIMDAddi8(x,y) _mm256_add_epi8(x,y)
#define SIMDAddi16(x,y) _mm256_add_epi16(x,y)
#define SIMDAddi32(x,y) _mm256_add_epi32(x,y)
#define SIMDAddi64(x,y) _mm256_add_epi64(x,y)
#define SIMDSub(x,y) _mm256_sub_ps(x,y)
#define SIMDSubi8(x,y) _mm256_sub_epi8(x,y)
#define SIMDSubi16(x,y) _mm256_sub_epi16(x,y)
#define SIMDSubi32(x,y) _mm256_sub_epi32(x,y)
#define SIMDSubi64(x,y) _mm256_sub_epi64(x,y)
#define SIMDMul(x,y) _mm256_mul_ps(x,y)
#define SIMDMuli(x,y) _mm256_mul_epi32(x,y)
#define SIMDAnd(x,y) _mm256_and_ps(x,y)
#define SIMDAndi(x,y) _mm256_and_si256(x,y)
#define SIMDAndNot(x,y) _mm256_andnot_ps(x,y)
#define SIMDAndNoti(x,y) _mm256_andnot_si256(x,y)
#define SIMDOr(x,y) _mm256_or_ps(x,y)
#define SIMDOri(x,y) _mm256_or_si256(x,y)
#define SIMDShiftLeft(a, n) (n) < 16 ? \
    _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), (16-(n))) : \
    _mm256_slli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), ((n)-16))

#define SIMDShiftRight(a, n) (n) < 16 ? \
    _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1)), (n)) : \
    _mm256_srli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1)), ((n)-16))
#define SIMDShiftLeftOnei16(x,y) _mm256_slli_epi16(x,y)
#define SIMDShiftLeftOnei32(x,y) _mm256_slli_epi32(x,y)
#define SIMDShiftLeftOnei64(x,y) _mm256_slli_epi64(x,y)
#define SIMDShiftRightOnei16(x,y) _mm256_srli_epi16(x,y)
#define SIMDShiftRightOnei32(x,y) _mm256_srli_epi32(x,y)
#define SIMDShiftRightOnei64(x,y) _mm256_srli_epi64(x,y)
#define SIMDEqual(x,y)  _mm256_cmp_ps(x,y,_CMP_EQ_OQ) 
#define SIMDEquali16(x,y) _mm256_cmpeq_epi16(x,y)
#define SIMDEquali8(x,y) _mm256_cmpeq_epi8(x,y)
#define SIMDEquali32(x,y) _mm256_cmpeq_epi32(x,y)
#define SIMDEquali64(x,y) _mm256_cmpeq_epi64(x,y)
#define SIMDGreaterThan(x,y) _mm256_cmp_ps(x,y,_CMP_GT_OQ)
#define SIMDGreaterThani16(x,y) _mm256_cmpgt_epi16(x,y)
#define SIMDGreaterThani8(x,y) _mm256_cmpgt_epi8(x,y)
#define SIMDGreaterThani32(x,y) _mm256_cmpgt_epi32(x,y)
#define SIMDGreaterThani64(x,y) _mm256_cmpgt_epi64(x,y) 
#define SIMDFloor(x) _mm256_floor_ps(x)
#define SIMDMax(x,y) _mm256_max_ps(x,y)
#define SIMDMaxi8(x,y) _mm256_max_epi8(x,y)
#define SIMDMaxi16(x,y) _mm256_max_epi16(x,y)
#define SIMDMaxi32(x,y) _mm256_max_epi32(x,y)
#define SIMDMaxi64(x,y) _mm256_max_epi64(x,y)
#define SIMDMin(x,y) _mm256_min_ps(x,y)
#define SIMDMini8(x,y) _mm256_min_epi8(x,y)
#define SIMDMini16(x,y) _mm256_min_epi16(x,y)
#define SIMDMini32(x,y) _mm256_min_epi32(x,y)

#define SIMDBlendV(x,y,z) _mm256_blendv_ps(x,y,z)
#define SIMDBlendVi8(x,y,z) _mm256_blendv_epi8(x,y,z)

// end of AVX2 only
 
#else // SSE4.1 SSE2

// start of SSE4.1 and SSE2
// m128 will be our base type
typedef __m128 SIMDf;   //for floats
typedef __m128i SIMDi; //for integers

#define SIMDStore(x,y) _mm_store_ps(x,y)
#define SIMDLoad(x) _mm_load_ps(x)
#define SIMDStorei(x,y) _mm_store_si128(x,y)
#define SIMDLoadi(x) _mm_load_si128(x)
#define SIMDSeti8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16) __mm_set_epi8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16)
#define SIMDSeti16(x1,x2,x3,x4,x5,x6,x7,x8) __mm_set_epi16(x1,x2,x3,x4,x5,x6,x7,x8)
#define SIMDSeti32(x1,x2,x3,x4) __mm_set_epi32(x1,x2,x3,x4)
#define SIMDSeti64(x,y) __mm_set_epi64(x,y)
#define SIMDSetOne(x) _mm_set1_ps(x)
#define SIMDSetZero() _mm_setzero_ps()
#define SIMDSetOnei8(x) _mm_set1_epi8(x)
#define SIMDSetOnei16(x) _mm_set1_epi16(x)
#define SIMDSetOnei32(x) _mm_set1_epi32(x)
#define SIMDSetOnei64(x) _mm_set1_epi64(x)
#define SIMDSetZeroi() _mm_setzero_si128()
#define SIMDAdd(x,y) _mm_add_ps(x,y)
#define SIMDAddi8(x,y) _mm_add_epi8(x,y)
#define SIMDAddi16(x,y) _mm_add_epi16(x,y)
#define SIMDAddi32(x,y) _mm_add_epi32(x,y)
#define SIMDAddi64(x,y) _mm_add_epi64(x,y)
#define SIMDSub(x,y) _mm_sub_ps(x,y)
#define SIMDSubi8(x,y) _mm_sub_epi8(x,y)
#define SIMDSubi16(x,y) _mm_sub_epi16(x,y)
#define SIMDSubi32(x,y) _mm_sub_epi32(x,y)
#define SIMDSubi64(x,y) _mm_sub_epi64(x,y)
#define SIMDMul(x,y) _mm_mul_ps(x,y)
#define SIMDMuli(x,y) _mm_mul_epi32(x,y)
#define SIMDAnd(x,y) _mm_and_ps(x,y)
#define SIMDAndi(x,y) _mm_and_si128(x,y)
#define SIMDAndNot(x,y) _mm_andnot_ps(x,y)
#define SIMDAndNoti(x,y) _mm_andnot_si128(x,y)
#define SIMDOr(x,y) _mm_or_ps(x,y)
#define SIMDOri(x,y) _mm_or_si128(x,y)
#define SIMDShiftLeft(x,y) _mm_slli_si128(x,y) // shift whole x by y bits
#define SIMDShiftRight(x,y) _mm_srli_si128(x,y)
#define SIMDShiftLeftOnei16(x,y) _mm_slli_epi16(x,y)
#define SIMDShiftLeftOnei32(x,y) _mm_slli_epi32(x,y)
#define SIMDShiftLeftOnei64(x,y) _mm_slli_epi64(x,y)
#define SIMDShiftRightOnei16(x,y) _mm_srli_epi16(x,y)
#define SIMDShiftRightOnei32(x,y) _mm_srli_epi32(x,y)
#define SIMDShiftRightOnei64(x,y) _mm_srli_epi64(x,y)
#define SIMDEqual(x,y)  _mm_cmpeq_ps(x,y)
#define SIMDEquali8(x,y) _mm_cmpeq_epi8(x,y)
#define SIMDEquali16(x,y) _mm_cmpeq_epi16(x,y)
#define SIMDEquali32(x,y) _mm_cmpeq_epi32(x,y)
#define SIMDGreaterThan(x,y) _mm_cmpgt_ps(x,y)
#define SIMDGreaterThani8(x,y) _mm_cmpgt_epi8(x,y)
#define SIMDGreaterThani16(x,y) _mm_cmpgt_epi16(x,y)
#define SIMDGreaterThani32(x,y) _mm_cmpgt_epi32(x,y)
#define SIMDLessThan(x,y) _mm_cmplt_ps(x,y)
#define SIMDLessThani8(x,y) _mm_cmplt_epi8(x,y) 
#define SIMDLessThani16(x,y) _mm_cmplt_epi16(x,y) 
#define SIMDLessThani32(x,y) _mm_cmplt_epi32(x,y) 
#define SIMDMax(x,y) _mm_max_ps(x,y)
#define SIMDMaxi16(x,y) _mm_max_epi16(x,y)
#define SIMDMin(x,y) _mm_min_ps(x,y)
#define SIMDMini16(x,y) _mm_min_epi16(x,y)

#define Maski16 __mmask8
#define Maski32 __mmask8

#ifdef __SSE4_1__

// start of SSE4.1 only
#define SIMDBlendV(x,y,z) _mm_blendv_ps(x,y,z)	    // z is __mask
#define SIMDBlendVi8(x,y,z) _mm_blendv_epi8(x,y,z)	
#define SIMDEquali64(x,y) _mm_cmpeq_epi64(x,y)
#define SIMDFloor(x) _mm_floor_ps(x)
#define SIMDMaxi8(x,y) _mm_max_epi8(x,y)
#define SIMDMini8(x,y) _mm_min_epi8(x,y)
#define SIMDMaxi32(x,y) _mm_max_epi32(x,y)
#define SIMDMini32(x,y) _mm_min_epi32(x,y)
// end of SSE4.1 only

#else  // SSE2

// start of SSE2 only
#define SIMDBlendV(x,y,z) SIMDOr(SIMDAndNot(z,x), SIMDAnd(z,y))   //if we don't have sse4
#define SIMDBlendVi8(x,y,z) SIMDOri(SIMDAndNoti(z,x), SIMDAndi(z,y))    //if we don't have sse4
#define SIMDMaxi8(x,y) SIMDBlendVi8(y, x, SIMDGreaterThani8(x,y))
#define SIMDMini8(x,y) SIMDBlendVi8(x, y, SIMDGreaterThani8(x,y))
#define SIMDMaxi32(x,y) SIMDBlendi32(y, x, SIMDGreaterThani32(x,y))
#define SIMDMini32(x,y) SIMDBlendi32(x, y, SIMDGreaterThani32((x,y))
// end of SSE2 only
// end of SSE4.1 and SSE2

#endif // SSE4.1

#endif // AVX2

// start of no AVX512BW (AVX2/SSE4.1/SSE2)
#define SIMDBlendi16(x,y,z) SIMDOri(SIMDAndNoti(z,x), SIMDAndi(z,y))
#define SIMDBlendi32(x,y,z) SIMDOri(SIMDAndNoti(z,x), SIMDAndi(z,y))
#define SIMDBlendi64(x,y,z) SIMDOri(SIMDAndNoti(z,x), SIMDAndi(z,y))

/* x = a == b ? c : d */ 
#define SIMDSetIfEquali8(x,a,b,c,d)  { x = SIMDBlendVi8(d, c, SIMDEquali8(a,b)); } 
#define SIMDSetIfEquali16(x,a,b,c,d) { x = SIMDBlendi16(d, c, SIMDEquali16(a,b)); } 
#define SIMDSetIfEquali32(x,a,b,c,d) { x = SIMDBlendi32(d, c, SIMDEquali32(a,b)); } 
#define SIMDSetIfEquali64(x,a,b,c,d) { x = SIMDBlendi64(d, c, SIMDEquali64(a,b)); } 
/* x = a > b ? c : d */
#define SIMDSetIfGreateri8(x,a,b,c,d)  { x = SIMDBlendVi8(d, c, SIMDGreaterThani8(a,b)); } 
#define SIMDSetIfGreateri16(x,a,b,c,d) { x = SIMDBlendi16(d, c, SIMDGreaterThani16(a,b)); } 
#define SIMDSetIfGreateri32(x,a,b,c,d) { x = SIMDBlendi32(d, c, SIMDGreaterThani32(a,b)); } 
#define SIMDSetIfGreateri64(x,a,b,c,d) { x = SIMDBlendi64(d, c, SIMDGreaterThani64(a,b)); } 
/* x = a < b ? c : d */
#define SIMDSetIfLessi8(x,a,b,c,d)  { x = BlendVi8(d, c, SIMDGreaterThani8(b,a)); } 
#define SIMDSetIfLessi16(x,a,b,c,d) { x = Blendi16(d, c, SIMDGreaterThani16(b,a)); } 
#define SIMDSetIfLessi32(x,a,b,c,d) { x = Blendi32(d, c, SIMDGreaterThani32(b,a)); } 
#define SIMDSetIfLessi64(x,a,b,c,d) { x = Blendi64(d, c, SIMDGreaterThani64(b,a)); } 

/* x = a > b ? c : d, y = a > b ? a : b */
#define SIMDGetIfGreateri8(x,y,a,b,c,d)  { SIMDi cmp = SIMDGreaterThani8(a,b);  x = SIMDBlendVi8(d, c, cmp); y = SIMDBlendVi8(b, a, cmp); } 
#define SIMDGetIfGreateri16(x,y,a,b,c,d) { SIMDi cmp = SIMDGreaterThani16(a,b); x = SIMDBlendi16(d, c, cmp); y = SIMDBlendi16(b, a, cmp); } 
#define SIMDGetIfGreateri32(x,y,a,b,c,d) { SIMDi cmp = SIMDGreaterThani32(a,b); x = SIMDBlendi32(d, c, cmp); y = SIMDBlendi32(b, a, cmp); } 
#define SIMDGetIfGreateri64(x,y,a,b,c,d) { SIMDi cmp = SIMDGreaterThani64(a,b); x = SIMDBlendi64(d, c, cmp); y = SIMDBlendi64(b, a, cmp); } 
/* x = a < b ? c : d, y = a < b ? a : b */
#define SIMDGetIfLessi8(x,y,a,b,c,d)  { SIMDi cmp = SIMDGreaterThani8(b,a);  x = SIMDBlendVi8(d, c, cmp); y = SIMDBlendVi8(b, a, cmp); } 
#define SIMDGetIfLessi16(x,y,a,b,c,d) { SIMDi cmp = SIMDGreaterThani16(b,a); x = SIMDBlendi16(d, c, cmp); y = SIMDBlendi16(b, a, cmp); } 
#define SIMDGetIfLessi32(x,y,a,b,c,d) { SIMDi cmp = SIMDGreaterThani32(b,a); x = SIMDBlendi32(d, c, cmp); y = SIMDBlendi32(b, a, cmp); } 
#define SIMDGetIfLessi64(x,y,a,b,c,d) { SIMDi cmp = SIMDGreaterThani64(b,a); x = SIMDBlendi64(d, c, cmp); y = SIMDBlendi64(b, a, cmp); } 
// end of no AVX512BW (AVX2/SSE4.1/SSE2)

#endif // AVX512BW

#ifdef __cplusplus
extern "C" {
#endif

// int simd_check(void);

/*
static void *SIMDMalloc(size_t size, size_t align) {
    void *ret = (void*)_mm_malloc(size, align);
    if (ret == NULL) {
        fprintf(stderr, "[%s] mm_Malloc fail!\nSize: %ld\n", __func__, size);
        exit(1);
    }
    else return ret;
}*/

// use posix_memalign
static void *SIMDMalloc(size_t size, size_t align) {
    void *ret; int res;
    res = posix_memalign(&ret, align, size);
    if (res != 0) {
        char error[10];
        if (res == EINVAL) strcpy(error, "EINVAR");
        else if (res == ENOMEM)
            strcpy(error, "ENOMEM");
        else strcpy(error, "Unknown");
        fprintf(stderr, "[%s] posix_memalign fail!\nSize: %ld, Error: %s\n", __func__, size, error);
        exit(1);
    }
    else return ret;
}
#ifdef __cplusplus
}
#endif


#endif // SIMD_INSTRUCTION_H
