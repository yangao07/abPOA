#ifndef SIMD_POA_ALIGN_H
#define SIMD_POA_ALIGN_H

int simd_bpoa_banded_global_align_sequence_with_graph(bpoa_graph_t *graph, uint8_t *query, int qlen, bpoa_para_t *bpt, int *n_cigar, bpoa_cigar_t **graph_cigar);

#endif
