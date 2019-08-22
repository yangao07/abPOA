#ifndef AGGLO_HIER_CLU_H
#define AGGLO_HIER_CLU_H

#ifdef __cplusplus
extern "C" {
#endif

int agglo_hier_clu(int **weight_matrix, int matrix_size, int clu_n, double min_fre, int **clu_ids, int *clu_ids_n);



#ifdef __cplusplus
}
#endif

#endif
