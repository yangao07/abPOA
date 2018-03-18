#ifndef SIMD_POA_ALIGN_H
#define SIMD_POA_ALIGN_H

int simd_abpoa_banded_global_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar);

#endif
