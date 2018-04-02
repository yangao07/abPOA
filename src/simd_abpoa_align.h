#ifndef SIMD_ABPOA_ALIGN_H
#define SIMD_ABPOA_ALIGN_H

int simd_abpoa_align_sequence_with_graph(abpoa_t *ab, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar);
abpoa_simd_matrix_t *abpoa_init_simd_matrix(void);
void abpoa_free_simd_matrix(abpoa_simd_matrix_t *abm);

#endif
