#ifndef ABPOA_ALIGN_H
#define ABPOA_ALIGN_H

#include "abpoa.h"
#include "abpoa_graph.h"

#define CHUNK_READ_N 1000

// SPOA
#define ABPOA_MATCH  2
#define ABPOA_MISMATCH  4
#define ABPOA_GAP_OPEN1  4
#define ABPOA_GAP_OPEN2  24
#define ABPOA_GAP_EXT1  2
#define ABPOA_GAP_EXT2  1

//#define ABPOA_MATCH  1
//#define ABPOA_MISMATCH  3
//#define ABPOA_GAP_OPEN  5
//#define ABPOA_GAP_EXT  2

#define ABPOA_GLOBAL_MODE 0
#define ABPOA_EXTEND_MODE 1
#define ABPOA_LOCAL_MODE 2
//#define ABPOA_SEMI_MODE 3

// gap mode
#define ABPOA_LINEAR_GAP 0
#define ABPOA_AFFINE_GAP 1
#define ABPOA_CONVEX_GAP 2

#define ABPOA_MULTIP     1
#define ABPOA_MIN_FRE    0.3

#define ABPOA_CIGAR_STR "MIDXSH"
#define ABPOA_CMATCH     0
#define ABPOA_CINS       1
#define ABPOA_CDEL       2
#define ABPOA_CDIFF      3
#define ABPOA_CSOFT_CLIP 4
#define ABPOA_CHARD_CLIP 5

#define HOP_OFF_SET   0
#define EOP_OFF_SET   2
#define FOP_OFF_SET   4

// calculate band range for each row:
// have: min_rank, max_rank, min_remain, max_remain
// then: min_len = min_rank + min_remain, max_len = max_rank + max_remain
// range: (min_of_two(min_rank, max_rank+qlen-max_len), max_of_two(min_rank+qlen-min_len, max_rank))
// with w: (min-w, max+w)
#define GET_DP_BEGIN(graph, w, i, qlen) MAX_OF_TWO(0,    MIN_OF_TWO(abpoa_graph_node_id_to_min_rank(graph, i), qlen - abpoa_graph_node_id_to_max_remain(graph, i))-w)
#define GET_DP_END(graph, w, i, qlen)   MIN_OF_TWO(qlen, MAX_OF_TWO(abpoa_graph_node_id_to_max_rank(graph, i), qlen - abpoa_graph_node_id_to_min_remain(graph, i))+w)

// calculate band range for each row:
// have: min_remain, max_remain, max_i (max in row)
// range: (min_of_two(max_i+1, qlen-max_remain), max_of_two(max_i+1, qlen-min_remain))
// with w: (min-w, max+w)
#define GET_AD_DP_BEGIN(graph, w, i, qlen) MAX_OF_TWO(0,    MIN_OF_TWO(abpoa_graph_node_id_to_min_rank(graph, i), qlen - abpoa_graph_node_id_to_max_remain(graph, i))-w)
#define GET_AD_DP_END(graph, w, i, qlen)   MIN_OF_TWO(qlen, MAX_OF_TWO(abpoa_graph_node_id_to_max_rank(graph, i), qlen - abpoa_graph_node_id_to_min_remain(graph, i))+w)

#define _set_max_score(best_score, best_i, best_j, score, i, j) { \
    if (score > best_score) {                                     \
        best_score = score; best_i = i; best_j = j;               \
    }                                                             \
}

#ifdef __cplusplus
extern "C" {
#endif

// TODO splice mode
/* banded global partial order graph alignment */
int abpoa_banded_global_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar);

/* Adaptive banded global partial order graph alignment */
int ada_abpoa_banded_global_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar);

/* Banded global partial order graph alignment */
int abpoa_global_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar);

static inline abpoa_cigar_t *abpoa_push_cigar(int *n_cigar, int *m_cigar, abpoa_cigar_t *cigar, int op, int len, int32_t node_id, int32_t query_id) {
    abpoa_cigar_t l = len;
    if (*n_cigar == 0 || (op != ABPOA_CINS && op != ABPOA_CSOFT_CLIP && op != ABPOA_CHARD_CLIP) || op != (cigar[(*n_cigar)-1] & 0xf)) {
        if (*n_cigar == *m_cigar) {
            *m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
            cigar = (abpoa_cigar_t*)_err_realloc(cigar, (*m_cigar) * sizeof(abpoa_cigar_t));
        }
        abpoa_cigar_t n_id = node_id, q_id = query_id;
        if (op == ABPOA_CMATCH || op == ABPOA_CDIFF) 
            cigar[(*n_cigar)++] = n_id << 34 | q_id << 4 | op;
        else if (op == ABPOA_CINS || op == ABPOA_CSOFT_CLIP || op == ABPOA_CHARD_CLIP) 
            cigar[(*n_cigar)++] = q_id << 34 | l << 4 | op;
        else if (op == ABPOA_CDEL)
            cigar[(*n_cigar)++] = n_id << 34 | l << 4 | op;
        else
            err_fatal(__func__, "Unknown cigar operation: %s\n", op);
    } else cigar[(*n_cigar)-1] += l << 4;

    return cigar;
}

static inline abpoa_cigar_t *abpoa_reverse_cigar(int n_cigar, abpoa_cigar_t *cigar) {
    int i; abpoa_cigar_t tmp;
    for (i = 0; i < n_cigar >> 1; ++i) {
        tmp = cigar[i];
        cigar[i] = cigar[n_cigar-1-i];
        cigar[n_cigar-1-i] = tmp;
    }
    return cigar;
}

static inline void abpoa_print_cigar(int n_cigar, abpoa_cigar_t *cigar, abpoa_graph_t *graph) {
    int i, op, len, node_id, query_id, index_i;
    int n[6] = {0, 0, 0, 0, 0, 0};
    for (i = 0; i < n_cigar; ++i) {
        op = cigar[i] & 0xf; node_id = (int)(cigar[i] >> 34); 
        len = query_id = (int)(cigar[i] >> 4) & 0x3fffffff;
        if (op == ABPOA_CMATCH || op == ABPOA_CDIFF) {
            index_i = abpoa_graph_node_id_to_index(graph, node_id);
            printf("1%c:%d,%d\t", ABPOA_CIGAR_STR[op], index_i, query_id);
            n[op] += 1;
        } else if (op == ABPOA_CDEL) {
            index_i = abpoa_graph_node_id_to_index(graph, node_id);
            printf("%d%c:%d\t", len, ABPOA_CIGAR_STR[op], index_i);
            n[op] += len;
        } else if (op == ABPOA_CINS || op == ABPOA_CSOFT_CLIP || op == ABPOA_CHARD_CLIP) { 
            printf("%d%c:%d\t", len, ABPOA_CIGAR_STR[op], node_id);
            n[op] += len;
        } else {
            err_fatal(__func__, "Unknown cigar operation: %s\n", op);
        }
    } printf("\n");
    for (i = 0; i < 6; ++i)
        printf("%d%c ", n[i], ABPOA_CIGAR_STR[i]);
    printf("\n");
}

// type of score matrix (DP_H, DP_E): 8/16/32
static inline void abpoa_backtrack(int *DP_H, int *DP_E, int matrix_col_n, int m,  int *mat, int8_t gap_e, int **pre_index, int *pre_n, uint8_t *backtrack_z, int best_i, int best_j, int z_col_n, abpoa_graph_t *graph, uint8_t *query, int *n_cigar, abpoa_cigar_t **graph_cigar) {
    int i, j, k, pre_i;
    if (n_cigar && graph_cigar) {
        int n_c = 0, s, m_c = 0, id, which, last_which;
        int op_shift[3] = {HOP_OFF_SET, EOP_OFF_SET, FOP_OFF_SET};
        uint8_t d;
        abpoa_cigar_t *cigar = 0;
        i = best_i, j = best_j, id = abpoa_graph_index_to_node_id(graph, i), last_which = 0;
        while (i > 0 && j > 0) {
            d = backtrack_z[(long)(i-1) * z_col_n + j-1];
            which = (d >> op_shift[last_which]) & 3;
            if (which == 0) { // match
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);
                s = mat[graph->node[id].base * m + query[j-1]];
                //s = graph->node[id].base == query[j-1] ? match : -mis;

                for (k = 0; k < pre_n[i]; ++k) {
                    pre_i = pre_index[i][k];
                    if (DP_H[pre_i * matrix_col_n + j-1] + s == DP_H[i * matrix_col_n + j]) {
                        i = pre_i;
                        break;
                    }
                }
                id = abpoa_graph_index_to_node_id(graph, i);
                j--; last_which = which;
            } else if (which == 1) { // deletion
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);
                if (last_which == 1) {
                    for (k = 0; k < pre_n[i]; ++k) {
                        pre_i = pre_index[i][k];
                        if (DP_E[pre_i * matrix_col_n + j] - gap_e == DP_E[i * matrix_col_n + j]) {
                            i = pre_i;
                            break;
                        }
                    }
                } else if (last_which == 0) {
                    for (k = 0; k < pre_n[i]; ++k) {
                        pre_i = pre_index[i][k];
                        if (DP_E[pre_i * matrix_col_n + j] == DP_H[i * matrix_col_n + j]) {
                            i = pre_i;
                            break;
                        }
                    }
                } else {
                    printf("\nunexpected cigar op.\n");
                }
                id = abpoa_graph_index_to_node_id(graph, i);
                last_which = which;
            } else { // insertion
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1);
                j--; last_which = which;
            }
        }
        if (j > 0) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, j, -1, j-1);
        // reverse cigar
        *graph_cigar = abpoa_reverse_cigar(n_c, cigar);
        *n_cigar = n_c;
#ifdef __DEBUG__
        abpoa_print_cigar(n_c, *graph_cigar, graph);
#endif
    }
}

#ifdef __cplusplus
}
#endif

#endif
