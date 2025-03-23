#ifdef ABPOA_SIMD_DISPATCH
#include <stdio.h>
#include <stdlib.h>
#include "abpoa.h"


#define SIMD_SSE     0x1
#define SIMD_SSE2    0x2
#define SIMD_SSE3    0x4
#define SIMD_SSSE3   0x8
#define SIMD_SSE4_1  0x10
#define SIMD_SSE4_2  0x20
#define SIMD_AVX     0x40
#define SIMD_AVX2    0x80
#define SIMD_AVX512F 0x100
#define SIMD_AVX512BW 0x200

#ifndef _MSC_VER
// adapted from https://github.com/01org/linux-sgx/blob/master/common/inc/internal/linux/cpuid_gnu.h
void __cpuidex(int cpuid[4], int func_id, int subfunc_id) {
#if defined(__x86_64__)
	__asm__ volatile ("cpuid"
			: "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#else // on 32bit, ebx can NOT be used as PIC code
	__asm__ volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
			: "=a" (cpuid[0]), "=r" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#endif
}
#endif

static int x86_simd(void) {
	int flag = 0, cpuid[4], max_id;
	__cpuidex(cpuid, 0, 0);
	max_id = cpuid[0];
	if (max_id == 0) return 0;
	__cpuidex(cpuid, 1, 0);
	if (cpuid[3]>>25&1) flag |= SIMD_SSE;
	if (cpuid[3]>>26&1) flag |= SIMD_SSE2;
	if (cpuid[2]>>0 &1) flag |= SIMD_SSE3;
	if (cpuid[2]>>9 &1) flag |= SIMD_SSSE3;
	if (cpuid[2]>>19&1) flag |= SIMD_SSE4_1;
	if (cpuid[2]>>20&1) flag |= SIMD_SSE4_2;
	if (cpuid[2]>>28&1) flag |= SIMD_AVX;
	if (max_id >= 7) {
		__cpuidex(cpuid, 7, 0);
		if (cpuid[1]>>5 &1) flag |= SIMD_AVX2;
		if (cpuid[1]>>16&1) flag |= SIMD_AVX512F;
		if (cpuid[1]>>30&1) flag |= SIMD_AVX512BW;
	}
	return flag;
}

static int simd_flag = -1;


int simd_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res) {
	extern void simd_sse2_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res);
	extern void simd_sse41_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res);
	extern void simd_avx2_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res);
	extern void simd_avx512_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res);
	if (simd_flag < 0) simd_flag = x86_simd();
	if (simd_flag & SIMD_AVX512BW) {
        // fprintf(stderr, "AVX512\n");
        simd_avx512_abpoa_align_sequence_to_subgraph(ab, abpt, beg_node_id, end_node_id, query, qlen, res);
    } else if (simd_flag & SIMD_AVX2) {
        // fprintf(stderr, "AVX2\n");
        simd_avx2_abpoa_align_sequence_to_subgraph(ab, abpt, beg_node_id, end_node_id, query, qlen, res);
    } else if (simd_flag & SIMD_SSE4_1) {
        // fprintf(stderr, "SSE41\n");
        simd_sse41_abpoa_align_sequence_to_subgraph(ab, abpt, beg_node_id, end_node_id, query, qlen, res);
    } else if (simd_flag & SIMD_SSE2) {
        // fprintf(stderr, "SSE2\n");
        simd_sse2_abpoa_align_sequence_to_subgraph(ab, abpt, beg_node_id, end_node_id, query, qlen, res);
    } else abort();
	return 0;
}
int simd_abpoa_align_sequence_to_graph(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, abpoa_res_t *res) {
	return simd_abpoa_align_sequence_to_subgraph(ab, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, query, qlen, res);
}

#endif
