#ifndef ABPOA_ALIGN_H
#define ABPOA_ALIGN_H

#include "abpoa.h"
#include "abpoa_graph.h"

#define CHUNK_READ_N 1024

#define ABPOA_MATCH  2
#define ABPOA_MISMATCH  4
#define ABPOA_GAP_OPEN1  4
#define ABPOA_GAP_OPEN2  24
#define ABPOA_GAP_EXT1  2
#define ABPOA_GAP_EXT2  1

#define ABPOA_M_OP   0x1
#define ABPOA_E1_OP  0x2
#define ABPOA_E2_OP  0x4
#define ABPOA_E_OP   0x6
#define ABPOA_F1_OP  0x8 
#define ABPOA_F2_OP  0x10
#define ABPOA_F_OP   0x18
#define ABPOA_ALL_OP 0x1f

#define DIPLOID_MIN_FREQ    0.3

// start and end of each band:
//   range: (min_of_two(max_left, qlen-remain), max_of_two(max_right, qlen-remain))
//   with extra band width: (range_min-w, range_max+w)
#define GET_AD_DP_BEGIN(graph, w, i, qlen) MAX_OF_TWO(0,    MIN_OF_TWO(abpoa_graph_node_id_to_max_pos_left(graph, i),  qlen - abpoa_graph_node_id_to_max_remain(graph, i))-w)
#define GET_AD_DP_END(graph, w, i, qlen)   MIN_OF_TWO(qlen, MAX_OF_TWO(abpoa_graph_node_id_to_max_pos_right(graph, i), qlen - abpoa_graph_node_id_to_max_remain(graph, i))+w)

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif
