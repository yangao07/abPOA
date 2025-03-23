#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "abpoa_align.h"
#include "abpoa_simd.h"
#include "simd_instruction.h"
#include "utils.h"

// #ifdef ABPOA_SIMD_DISPATCH
// #if __AVX512BW__
// SIMD_para_t _simd_p8_avx512  = {512,  8, 6, 64, 64, -1};
// SIMD_para_t _simd_p16_avx512 = {512, 16, 5, 32, 64, -1};
// SIMD_para_t _simd_p32_avx512 = {512, 32, 4, 16, 64, -1};
// SIMD_para_t _simd_p64_avx512 = {512, 64, 3,  8, 64, -1};
// #elif defined(__AVX2__)
// SIMD_para_t _simd_p8_avx2  = {256,  8, 5, 32, 32, -1};
// SIMD_para_t _simd_p16_avx2 = {256, 16, 4, 16, 32, -1};
// SIMD_para_t _simd_p32_avx2 = {256, 32, 3,  8, 32, -1};
// SIMD_para_t _simd_p64_avx2 = {256, 64, 2,  4, 32, -1};
// #else
// SIMD_para_t _simd_p8_sse  = {128,  8, 4, 16, 16, -1};
// SIMD_para_t _simd_p16_sse = {128, 16, 3,  8, 16, -1};
// SIMD_para_t _simd_p32_sse = {128, 32, 2,  4, 16, -1};
// SIMD_para_t _simd_p64_sse = {128, 64, 1,  2, 16, -1};
// #endif
// #else
// #if __AVX512BW__
// SIMD_para_t _simd_p8  = {512,  8, 6, 64, 64, -1};
// SIMD_para_t _simd_p16 = {512, 16, 5, 32, 64, -1};
// SIMD_para_t _simd_p32 = {512, 32, 4, 16, 64, -1};
// SIMD_para_t _simd_p64 = {512, 64, 3,  8, 64, -1};
// #elif defined(__AVX2__)
// SIMD_para_t _simd_p8  = {256,  8, 5, 32, 32, -1};
// SIMD_para_t _simd_p16 = {256, 16, 4, 16, 32, -1};
// SIMD_para_t _simd_p32 = {256, 32, 3,  8, 32, -1};
// SIMD_para_t _simd_p64 = {256, 64, 2,  4, 32, -1};
// #else
// SIMD_para_t _simd_p8  = {128,  8, 4, 16, 16, -1};
// SIMD_para_t _simd_p16 = {128, 16, 3,  8, 16, -1};
// SIMD_para_t _simd_p32 = {128, 32, 2,  4, 16, -1};
// SIMD_para_t _simd_p64 = {128, 64, 1,  2, 16, -1};
// #endif
// #endif


#define print_simd(s, str, score_t) {                                \
    int _i; score_t *_a = (score_t*)(s);                             \
    fprintf(stderr, "%s\t", str);                                    \
    for (_i = 0; _i < SIMDTotalBytes / (int)sizeof(score_t); ++_i) { \
        fprintf(stderr, "%d\t", _a[_i]);                             \
    } fprintf(stderr, "\n");                                         \
}

#define simd_abpoa_print_lg_matrix(score_t, beg_index, end_index) { \
    for (j = 0; j < end_index-beg_index; ++j) {                     \
        fprintf(stderr, "index: %d\t", j);                          \
        dp_h = DP_H + j * dp_sn;                                    \
        _dp_h = (score_t*)dp_h;                                     \
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {                  \
            fprintf(stderr, "%d:(%d)\t", i, _dp_h[i]);              \
        } fprintf(stderr, "\n");                                    \
    }                                                               \
}

#define simd_abpoa_print_ag_matrix(score_t, beg_index, end_index) {  \
    for (j = 0; j < end_index-beg_index; ++j) {                      \
        fprintf(stderr, "index: %d\t", j);                           \
        dp_h = DP_HEF + j * 3 * dp_sn; dp_e1 = dp_h + dp_sn;         \
        _dp_h = (score_t*)dp_h, _dp_e1 = (score_t*)dp_e1;            \
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {                   \
            fprintf(stderr, "%d:(%d,%d)\t", i, _dp_h[i], _dp_e1[i]); \
        } fprintf(stderr, "\n");                                     \
    }                                                                \
}

#define debug_simd_abpoa_print_cg_matrix_row(str, score_t, index_i) {                                    \
    score_t *_dp_h = (score_t*)dp_h, *_dp_e1 = (score_t*)dp_e1;                                          \
    score_t *_dp_e2 = (score_t*)dp_e2, *_dp_f1 = (score_t*)dp_f1, *_dp_f2 = (score_t*)dp_f2;             \
    fprintf(stderr, "%s\tindex: %d\t", str, index_i);                                                    \
    for (i = dp_beg[index_i]; i <= dp_end[index_i]; ++i) {                                               \
        fprintf(stderr, "%d:(%d,%d,%d,%d,%d)\t", i, _dp_h[i], _dp_e1[i],_dp_e2[i], _dp_f1[i],_dp_f2[i]); \
    } fprintf(stderr, "\n");                                                                             \
}

#define simd_abpoa_print_cg_matrix(score_t, beg_index, end_index) {                                          \
    for (j = 0; j < end_index-beg_index; ++j) {                                                              \
        fprintf(stderr, "index: %d\t", j);                                                                   \
        dp_h=DP_H2E2F+j*5*dp_sn; dp_e1=dp_h+dp_sn; dp_e2=dp_e1+dp_sn; dp_f1=dp_e2+dp_sn; dp_f2=dp_f1+dp_sn;  \
        score_t *_dp_h=(score_t*)dp_h, *_dp_e1=(score_t*)dp_e1, *_dp_e2=(score_t*)dp_e2;                     \
        score_t *_dp_f1=(score_t*)dp_f1, *_dp_f2=(score_t*)dp_f2;                                            \
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {                                                           \
            fprintf(stderr, "%d:(%d,%d,%d,%d,%d)\t", i, _dp_h[i], _dp_e1[i],_dp_e2[i], _dp_f1[i],_dp_f2[i]); \
        } fprintf(stderr, "\n");                                                                             \
    }                                                                                                        \
}

/* max_pos_left/right: left/right boundary of max column index for each row, based on the pre_nodes' DP score */
/* === workflow of alignment === */
/* a. global:
 * (1) alloc mem
 * (2) init for first row
 * (3) DP for each row
 * (3.2) if use_ada, update max_pos_left/right
 * (4) find best_i/j, backtrack
 * b. extend:
 * (1) alloc mem
 * (2) init for first row
 * (3) DP for each row
 * (3.2) find max of current row
 * (3.3) z-drop, set_max_score
 * (3.4) if use_ada, update max_pos_left/right
 */

// backtrack order:
// Match/Mismatch, Deletion, Insertion
#define simd_abpoa_lg_backtrack(score_t) {                                                                  \
    int i, j, k, pre_i, n_c = 0, s, is_match, m_c = 0, hit, id, _start_i, _start_j;                         \
    SIMDi *dp_h; score_t *_dp_h=NULL, *_pre_dp_h; abpoa_cigar_t *cigar = 0;                                 \
    i = best_i, j = best_j, _start_i = best_i, _start_j = best_j;                                           \
    id = abpoa_graph_index_to_node_id(graph, i+beg_index);                                                  \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, qlen-j, -1, qlen-1);         \
    dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h; int path_score = 0;                                    \
    int look_for_first_gap_at_end = abpt->put_gap_at_end; /* prefer to keep the last gap at end (begining of the backtrack) */ \
    int put_gap_on_right = abpt->put_gap_on_right;                                                          \
    while (i > 0 && j > 0) {                                                                                \
        if (abpt->align_mode == ABPOA_LOCAL_MODE && _dp_h[j] == 0) break;                                   \
        _start_i = i, _start_j = j;                                                                         \
        int *pre_index_i = pre_index[i];                                                                    \
        s = mat[m * graph->node[id].base + query[j-1]]; hit = 0;                                            \
        is_match = graph->node[id].base == query[j-1];                                                      \
        if (put_gap_on_right == 0 && look_for_first_gap_at_end == 0) {                                      \
            for (k = 0; k < pre_n[i]; ++k) { /* match/mismatch */                                           \
                pre_i = pre_index_i[k];                                                                     \
                if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);            \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_H + pre_i * dp_sn);                                               \
                if (_pre_dp_h[j-1] + s + path_score == _dp_h[j]) {                                          \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; hit = 1; id = abpoa_graph_index_to_node_id(graph, i+beg_index);         \
                    dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h;                                        \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) { /* deletion */                                                                      \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);            \
                if (j < dp_beg[pre_i] || j > dp_end[pre_i]) continue;                                       \
                _pre_dp_h = (score_t*)( DP_H + pre_i * dp_sn);                                              \
                if (_pre_dp_h[j] - gap_ext1 + path_score == _dp_h[j]) {                                     \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                    \
                    i = pre_i; hit = 1; id = abpoa_graph_index_to_node_id(graph, i+beg_index);              \
                    dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h;                                        \
                    if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                           \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) { /* insertion */                                                                     \
            if (_dp_h[j-1] - gap_ext1 == _dp_h[j]) {                                                        \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1); j--;                   \
                if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                               \
                hit = 1; ++res->n_aln_bases;                                                                \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) { /* match/mismatch */                                                                \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);            \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_H + pre_i * dp_sn);                                               \
                if (_pre_dp_h[j-1] + s + path_score == _dp_h[j]) { /* match/mismatch */                     \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; hit = 1; id = abpoa_graph_index_to_node_id(graph, i+beg_index);         \
                    dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h;                                        \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    look_for_first_gap_at_end = 0;                                                          \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) err_fatal_simple("Error in lg_backtrack.");                                           \
        simd_output_pre_nodes(pre_index[i], pre_n[i], i, j, 0, abpt->verbose);                              \
    }                                                                                                       \
    if (j > 0) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, j, -1, j-1);                         \
    /* reverse cigar */                                                                                     \
    res->graph_cigar = abpt->rev_cigar ? cigar : abpoa_reverse_cigar(n_c, cigar);                           \
    res->n_cigar = n_c; res->m_cigar = m_c;                                                                 \
    res->node_e = abpoa_graph_index_to_node_id(graph, best_i+beg_index), res->query_e=best_j-1; /*0-based*/ \
    res->node_s = abpoa_graph_index_to_node_id(graph, _start_i+beg_index), res->query_s=_start_j-1;         \
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/                                                        \
}

#define simd_abpoa_ag_backtrack(score_t) {                                                                  \
    int i, j, k, pre_i, n_c = 0, s, is_match, m_c = 0, id, hit, cur_op = ABPOA_ALL_OP, _start_i, _start_j;  \
    score_t *_dp_h, *_dp_e1, *_dp_f1, *_pre_dp_h, *_pre_dp_e1; abpoa_cigar_t *cigar = 0;                    \
    i = best_i, j = best_j; _start_i = best_i, _start_j = best_j;                                           \
    id = abpoa_graph_index_to_node_id(graph, i+beg_index);                                                  \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, qlen-j, -1, qlen-1);         \
    SIMDi *dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h; int path_score = 0;                       \
    int look_for_first_gap_at_end = abpt->put_gap_at_end; /* prefer to keep the last gap at end (begining of the backtrack) */ \
    int put_gap_on_right = abpt->put_gap_on_right;                                                          \
    while (i > 0 && j > 0) {                                                                                \
        if (abpt->align_mode == ABPOA_LOCAL_MODE && _dp_h[j] == 0) break;                                   \
        _start_i = i, _start_j = j;                                                                         \
        int *pre_index_i = pre_index[i];                                                                    \
        s = mat[m * graph->node[id].base + query[j-1]]; hit = 0;                                            \
        is_match = graph->node[id].base == query[j-1];                                                      \
        if (put_gap_on_right == 0 && look_for_first_gap_at_end == 0) {                                      \
            if (cur_op & ABPOA_M_OP) { /* match/mismatch */                                                 \
                for (k = 0; k < pre_n[i]; ++k) {                                                            \
                    pre_i = pre_index_i[k];                                                                 \
                    if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);        \
                    if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                               \
                    _pre_dp_h = (score_t*)(DP_HEF + dp_sn * pre_i * 3);                                     \
                    if (_pre_dp_h[j-1] + s + path_score == _dp_h[j]) {                                      \
                        cur_op = ABPOA_ALL_OP; hit = 1;                                                     \
                        cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);              \
                        i = pre_i; --j; id = abpoa_graph_index_to_node_id(graph, i+beg_index);              \
                        dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                              \
                        ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                       \
                        break;                                                                              \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_E1_OP) { /* deletion */                                              \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);            \
                if (j < dp_beg[pre_i] || j > dp_end[pre_i]) continue;                                       \
                _pre_dp_e1 = (score_t*)(DP_HEF + dp_sn * (pre_i * 3 + 1));                                  \
                if (cur_op & ABPOA_M_OP) {                                                                  \
                    if (_dp_h[j] == _pre_dp_e1[j] + path_score) {                                           \
                        _pre_dp_h = (score_t*)(DP_HEF + dp_sn * pre_i * 3);                                 \
                        if (_pre_dp_h[j] - gap_oe1 == _pre_dp_e1[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;      \
                        else cur_op = ABPOA_E1_OP;                                                          \
                        hit = 1;                                                                            \
                        cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                \
                        i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);                   \
                        dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                              \
                        if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                       \
                        break;                                                                              \
                    }                                                                                       \
                } else {                                                                                    \
                    _dp_e1 = (score_t*)(dp_h + dp_sn);                                                      \
                    if (_dp_e1[j] == _pre_dp_e1[j] - gap_ext1 + path_score) {                               \
                        _pre_dp_h = (score_t*)(DP_HEF + dp_sn * pre_i * 3);                                 \
                        if (_pre_dp_h[j] - gap_oe1 == _pre_dp_e1[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;      \
                        else cur_op = ABPOA_E1_OP;                                                          \
                        hit = 1;                                                                            \
                        cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                \
                        i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);                   \
                        dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                              \
                        if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                       \
                        break;                                                                              \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_F_OP) { /* insertion */                                              \
            _dp_f1 = (score_t*)(dp_h + dp_sn * 2);                                                          \
            if (cur_op & ABPOA_M_OP) {                                                                      \
                if (_dp_h[j] == _dp_f1[j]) {                                                                \
                    if (_dp_h[j-1] - gap_oe1 == _dp_f1[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;       \
                    else if (_dp_f1[j-1] - gap_ext1 == _dp_f1[j]) cur_op = ABPOA_F1_OP, hit = 1;            \
                }                                                                                           \
            } else {                                                                                        \
                if (_dp_h[j-1] - gap_oe1 == _dp_f1[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;           \
                else if (_dp_f1[j-1] - gap_ext1 == _dp_f1[j]) cur_op = ABPOA_F1_OP, hit = 1;                \
            }                                                                                               \
            if (hit == 1) {                                                                                 \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1); --j;                   \
                if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                               \
                ++res->n_aln_bases;                                                                         \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_M_OP) {                                                              \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);            \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_HEF + dp_sn * pre_i * 3);                                         \
                if (_pre_dp_h[j-1] + s + path_score == _dp_h[j]) {                                          \
                    cur_op = ABPOA_ALL_OP; hit = 1;                                                         \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; id = abpoa_graph_index_to_node_id(graph, i+beg_index);                  \
                    dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                                  \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    look_for_first_gap_at_end = 0;                                                          \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) err_fatal_simple("Error in ag_backtrack.");                                           \
        simd_output_pre_nodes(pre_index[i], pre_n[i], i, j, cur_op, abpt->verbose);                         \
    }                                                                                                       \
    if (j > 0) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, j, -1, j-1);                         \
    /* reverse cigar */                                                                                     \
    res->graph_cigar = abpt->rev_cigar ? cigar : abpoa_reverse_cigar(n_c, cigar);                           \
    res->n_cigar = n_c; res->m_cigar = m_c;                                                                 \
    res->node_e = abpoa_graph_index_to_node_id(graph, best_i+beg_index), res->query_e=best_j-1; /*0-based*/ \
    res->node_s = abpoa_graph_index_to_node_id(graph, _start_i+beg_index), res->query_s=_start_j-1;         \
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/                                                        \
}

#define simd_abpoa_cg_backtrack(score_t) {                                                                  \
    int i, j, k, pre_i, n_c = 0, s, is_match, m_c = 0, id, hit, cur_op = ABPOA_ALL_OP, _start_i, _start_j;  \
    score_t *_dp_h, *_dp_e1, *_dp_e2, *_dp_f1, *_dp_f2, *_pre_dp_h, *_pre_dp_e1, *_pre_dp_e2;               \
    abpoa_cigar_t *cigar = 0;                                                                               \
    i = best_i, j = best_j, _start_i = best_i, _start_j = best_j;                                           \
    id = abpoa_graph_index_to_node_id(graph, i+beg_index);                                                  \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, qlen-best_j, -1, qlen-1);    \
    SIMDi *dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h; int path_score = 0;                     \
    int look_for_first_gap_at_end = abpt->put_gap_at_end; /* prefer to keep the last gap at the end (begining of backtrack) */ \
    int put_gap_on_right = abpt->put_gap_on_right; /* prefer to keep gaps at the left-most position, minimap2-like */ \
    while (i > 0 && j > 0) {                                                                                \
        if (abpt->align_mode == ABPOA_LOCAL_MODE && _dp_h[j] == 0) break;                                   \
        _start_i = i, _start_j = j;                                                                         \
        int *pre_index_i = pre_index[i];                                                                    \
        s = mat[m * graph->node[id].base + query[j-1]]; hit = 0;                                            \
        is_match = graph->node[id].base == query[j-1];                                                      \
        if (put_gap_on_right == 0 && look_for_first_gap_at_end == 0) {                                      \
            if (cur_op & ABPOA_M_OP) { /* match/mismatch */                                                 \
                for (k = 0; k < pre_n[i]; ++k) {                                                            \
                    pre_i = pre_index_i[k];                                                                 \
                    if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);        \
                    if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                               \
                    _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                                   \
                    if (_pre_dp_h[j-1] + s + path_score == _dp_h[j]) {                                      \
                        cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);              \
                        i = pre_i; --j; id = abpoa_graph_index_to_node_id(graph, i+beg_index); hit = 1;     \
                        dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                            \
                        cur_op = ABPOA_ALL_OP;                                                              \
                        ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                       \
                        break;                                                                              \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_E_OP) { /* deletion */                                               \
            _dp_e1 = (score_t*)(dp_h + dp_sn); _dp_e2 = (score_t*)(dp_h + dp_sn * 2);                       \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);            \
                if (j < dp_beg[pre_i] || j > dp_end[pre_i]) continue;                                       \
                _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                                       \
                if (cur_op & ABPOA_E1_OP) {                                                                 \
                    _pre_dp_e1 = (score_t*)(DP_H2E2F + dp_sn * (pre_i * 5 + 1));                            \
                    if (cur_op & ABPOA_M_OP) {                                                              \
                        if (_dp_h[j] == _pre_dp_e1[j] + path_score) {                                       \
                            if (_pre_dp_h[j] - gap_oe1 == _pre_dp_e1[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;  \
                            else cur_op = ABPOA_E1_OP;                                                      \
                            hit = 1; cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);   \
                            i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);               \
                            dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                        \
                            if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                   \
                            break;                                                                          \
                        }                                                                                   \
                    } else {                                                                                \
                        if (_dp_e1[j] == _pre_dp_e1[j] - gap_ext1 + path_score) {                           \
                            if (_pre_dp_h[j] - gap_oe1 == _pre_dp_e1[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;  \
                            else cur_op = ABPOA_E1_OP;                                                      \
                            hit = 1; cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);   \
                            i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);               \
                            dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                        \
                            if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                   \
                            break;                                                                          \
                        }                                                                                   \
                    }                                                                                       \
                }                                                                                           \
                if (cur_op & ABPOA_E2_OP) {                                                                 \
                    _pre_dp_e2 = (score_t*)(DP_H2E2F + dp_sn * (pre_i * 5 + 2));                            \
                    if (cur_op & ABPOA_M_OP) {                                                              \
                        if (_dp_h[j] == _pre_dp_e2[j] + path_score) {                                       \
                            if (_pre_dp_h[j] - gap_oe2 == _pre_dp_e2[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;  \
                            else cur_op = ABPOA_E2_OP;                                                      \
                            hit = 1; cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);   \
                            i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);               \
                            dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                        \
                            if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                   \
                            break;                                                                          \
                        }                                                                                   \
                    } else {                                                                                \
                        if (_dp_e2[j] == _pre_dp_e2[j] - gap_ext2 + path_score) {                           \
                            if (_pre_dp_h[j] - gap_oe2 == _pre_dp_e2[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;  \
                            else cur_op = ABPOA_E2_OP;                                                      \
                            hit = 1; cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);   \
                            i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);               \
                            dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                        \
                            if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                   \
                            break;                                                                          \
                        }                                                                                   \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_F_OP) { /* insertion */                                              \
            if (cur_op & ABPOA_F1_OP) {                                                                     \
                _dp_f1 = (score_t*)(dp_h + dp_sn * 3);                                                      \
                if (cur_op & ABPOA_M_OP) {                                                                  \
                    if (_dp_h[j] == _dp_f1[j]) {                                                            \
                        if (_dp_h[j-1] - gap_oe1 == _dp_f1[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;   \
                        else if (_dp_f1[j-1] - gap_ext1 == _dp_f1[j]) cur_op = ABPOA_F1_OP, hit = 1;        \
                    }                                                                                       \
                } else {                                                                                    \
                    if (_dp_h[j-1] - gap_oe1 == _dp_f1[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;       \
                    else if (_dp_f1[j-1] - gap_ext1 == _dp_f1[j]) cur_op = ABPOA_F1_OP, hit = 1;            \
                }                                                                                           \
            }                                                                                               \
            if (hit == 0 && cur_op & ABPOA_F2_OP) {                                                         \
                _dp_f2 = (score_t*)(dp_h + dp_sn * 4);                                                      \
                if (cur_op & ABPOA_M_OP) {                                                                  \
                    if (_dp_h[j] == _dp_f2[j]) {                                                            \
                        if (_dp_h[j-1] - gap_oe2 == _dp_f2[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;   \
                        else if (_dp_f2[j-1] - gap_ext2 == _dp_f2[j]) cur_op = ABPOA_F2_OP, hit = 1;        \
                    }                                                                                       \
                } else {                                                                                    \
                    if (_dp_h[j-1] - gap_oe2 == _dp_f2[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;       \
                    else if (_dp_f2[j-1] - gap_ext2 == _dp_f2[j]) cur_op = ABPOA_F2_OP, hit = 1;            \
                }                                                                                           \
            }                                                                                               \
            if (hit == 1) {                                                                                 \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1); --j;                   \
                if (look_for_first_gap_at_end) look_for_first_gap_at_end = 0;                               \
                ++res->n_aln_bases;                                                                         \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && (cur_op & ABPOA_M_OP)) { /* match/mismatch */                                       \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, id, k);            \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                                       \
                if (_pre_dp_h[j-1] + s + path_score == _dp_h[j]) {                                          \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; id = abpoa_graph_index_to_node_id(graph, i+beg_index); hit = 1;         \
                    dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                                \
                    cur_op = ABPOA_ALL_OP;                                                                  \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    look_for_first_gap_at_end = 0;                                                          \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) err_fatal_simple("Error in cg_backtrack.");                                           \
        simd_output_pre_nodes(pre_index[i], pre_n[i], i, j, cur_op, abpt->verbose);                         \
    }                                                                                                       \
    if (j > 0) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, j, -1, j-1);                         \
    /* reverse cigar */                                                                                     \
    res->graph_cigar = abpt->rev_cigar ? cigar : abpoa_reverse_cigar(n_c, cigar);                           \
    res->n_cigar = n_c; res->m_cigar = m_c;                                                                 \
    res->node_e = abpoa_graph_index_to_node_id(graph, best_i+beg_index), res->query_e=best_j-1; /*0-based*/ \
    res->node_s = abpoa_graph_index_to_node_id(graph, _start_i+beg_index), res->query_s=_start_j-1;         \
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/                                                        \
}

// simd_abpoa_va
// simd_abpoa_ag_only_var
// sim_abpoa_init_var
#define simd_abpoa_var(score_t, sp, SIMDSetOne, SIMDShiftOneN)                                      \
    /* int tot_dp_sn = 0; */                                                                        \
    abpoa_graph_t *graph = ab->abg; abpoa_simd_matrix_t *abm = ab->abm;                             \
    int64_t matrix_row_n = end_index-beg_index+1, matrix_col_n = qlen + 1;                          \
    int **pre_index, *pre_n, _pre_index, _pre_n, pre_i;                                             \
    int i, j, k, *dp_beg, *dp_beg_sn, *dp_end, *dp_end_sn, node_id, index_i, dp_i;                  \
    int beg, end, beg_sn, end_sn, _beg_sn, _end_sn, pre_beg_sn, pre_end, sn_i;                      \
    int pn, log_n, size; int64_t qp_sn, dp_sn; /* pn: # value per SIMDi, qp_sn/dp_sn: segmented length*/ \
    SIMDi *dp_h, *pre_dp_h, *qp, *qi=NULL;                                                          \
    score_t *_dp_h=NULL, *_qi, best_score = sp.inf_min, inf_min = sp.inf_min;                       \
    int *mat = abpt->mat, m = abpt->m; score_t gap_ext1 = abpt->gap_ext1;                           \
    int w = abpt->wb < 0 ? qlen : abpt->wb+(int)(abpt->wf*qlen); /* when w < 0, do whole global */  \
    int best_i = 0, best_j = 0, best_id = 0, max, left_max_i=-1, right_max_i=-1;                    \
    SIMDi zero = SIMDSetZeroi(), SIMD_INF_MIN = SIMDSetOne(inf_min);                                \
    pn = sp.num_of_value; qp_sn = dp_sn = (matrix_col_n + pn - 1) / pn;                             \
    log_n = sp.log_num, size = sp.size; qp = abm->s_mem;                                            \
    int set_num; SIMDi *PRE_MASK, *SUF_MIN, *PRE_MIN;                                               \
    PRE_MASK = (SIMDi*)SIMDMalloc((pn+1) * size, size);                                             \
    SUF_MIN = (SIMDi*)SIMDMalloc((pn+1) * size, size);                                              \
    PRE_MIN = (SIMDi*)SIMDMalloc(pn * size, size);                                                  \
    for (i = 0; i < pn; ++i) {                                                                      \
        score_t *pre_mask = (score_t*)(PRE_MASK+i);                                                 \
        for (j = 0; j <= i; ++j) pre_mask[j] = -1;                                                  \
        for (j = i+1; j < pn; ++j) pre_mask[j] = 0;                                                 \
    } PRE_MASK[pn] = PRE_MASK[pn-1];                                                                \
    SUF_MIN[0] = SIMDShiftLeft(SIMD_INF_MIN, SIMDShiftOneN);                                        \
    for (i = 1; i < pn; ++i)                                                                        \
        SUF_MIN[i] = SIMDShiftLeft(SUF_MIN[i-1], SIMDShiftOneN); SUF_MIN[pn] = SUF_MIN[pn-1];       \
    for (i = 1; i < pn; ++i) {                                                                      \
        score_t *pre_min = (score_t*)(PRE_MIN + i);                                                 \
        for (j = 0; j < i; ++j) pre_min[j] = inf_min;                                               \
        for (j = i; j < pn; ++j) pre_min[j] = 0;                                                    \
    }

#define simd_abpoa_lg_only_var(score_t, SIMDSetOne, SIMDAdd)              \
    SIMDi *DP_H = qp + qp_sn * abpt->m; qi = DP_H + dp_sn * matrix_row_n; \
    SIMDi GAP_E1 = SIMDSetOne(gap_ext1);                                  \
    SIMDi *GAP_E1S =  (SIMDi*)SIMDMalloc(log_n * size, size);             \
    GAP_E1S[0] = GAP_E1;                                                  \
    for (i = 1; i < log_n; ++i) {                                         \
        GAP_E1S[i] = SIMDAdd(GAP_E1S[i-1], GAP_E1S[i-1]);                 \
    }

#define simd_abpoa_ag_only_var(score_t, SIMDSetOne, SIMDAdd)                                            \
    score_t *_dp_e1, *_dp_f1, gap_open1 = abpt->gap_open1, gap_oe1 = abpt->gap_open1 + abpt->gap_ext1;  \
    SIMDi *DP_HEF, *dp_e1, *pre_dp_e1, *dp_f1; int pre_end_sn;                                          \
    DP_HEF = qp + qp_sn * abpt->m; qi = DP_HEF + dp_sn * matrix_row_n * 3;                              \
    SIMDi GAP_O1 = SIMDSetOne(gap_open1), GAP_E1 = SIMDSetOne(gap_ext1), GAP_OE1 = SIMDSetOne(gap_oe1); \
    SIMDi *GAP_E1S =  (SIMDi*)SIMDMalloc(log_n * size, size);  GAP_E1S[0] = GAP_E1;                     \
    for (i = 1; i < log_n; ++i) {                                                                       \
        GAP_E1S[i] = SIMDAdd(GAP_E1S[i-1], GAP_E1S[i-1]);                                               \
    }

#define simd_abpoa_cg_only_var(score_t, SIMDSetOne, SIMDAdd)                                                      \
    score_t *_dp_e1, *_dp_e2, *_dp_f1, *_dp_f2, gap_open1 = abpt->gap_open1, gap_oe1 = gap_open1 + gap_ext1;      \
    score_t gap_open2 = abpt->gap_open2, gap_ext2 = abpt->gap_ext2, gap_oe2 = gap_open2 + gap_ext2;               \
    SIMDi *DP_H2E2F, *dp_e1, *dp_e2, *dp_f1, *dp_f2, *pre_dp_e1, *pre_dp_e2; int pre_end_sn;                      \
    SIMDi GAP_O1 = SIMDSetOne(gap_open1), GAP_O2 = SIMDSetOne(gap_open2);                                         \
    SIMDi GAP_E1 = SIMDSetOne(gap_ext1), GAP_E2 = SIMDSetOne(gap_ext2);                                           \
    SIMDi GAP_OE1 = SIMDSetOne(gap_oe1), GAP_OE2 = SIMDSetOne(gap_oe2);                                           \
    DP_H2E2F = qp + qp_sn * abpt->m; qi = DP_H2E2F + dp_sn * matrix_row_n * 5;                                    \
    SIMDi *GAP_E1S =  (SIMDi*)SIMDMalloc(log_n * size, size), *GAP_E2S =  (SIMDi*)SIMDMalloc(log_n * size, size); \
    GAP_E1S[0] = GAP_E1; GAP_E2S[0] = GAP_E2;                                                                     \
    for (i = 1; i < log_n; ++i) {                                                                                 \
        GAP_E1S[i] = SIMDAdd(GAP_E1S[i-1], GAP_E1S[i-1]);                                                         \
        GAP_E2S[i] = SIMDAdd(GAP_E2S[i-1], GAP_E2S[i-1]);                                                         \
    }

#define simd_abpoa_init_var(score_t) {                                                             \
    /* generate the query profile */                                                               \
    for (i = 0; i < qp_sn * abpt->m; ++i) qp[i] = SIMD_INF_MIN;                                    \
    for (k = 0; k < abpt->m; ++k) { /* SIMD parallelization */                                     \
        int *p = &mat[k * abpt->m];                                                                \
        score_t *_qp = (score_t*)(qp + k * qp_sn); _qp[0] = 0;                                     \
        for (j = 0; j < qlen; ++j) _qp[j+1] = (score_t)p[query[j]];                                \
        for (j = qlen+1; j < qp_sn * pn; ++j) _qp[j] = 0;                                          \
    }                                                                                              \
    if (abpt->wb>=0 || abpt->align_mode==ABPOA_LOCAL_MODE || abpt->align_mode==ABPOA_EXTEND_MODE){ \
        _qi = (score_t*)qi; /* query index */                                                      \
        for (i = 0; i <= qlen; ++i) _qi[i] = i;                                                    \
        for (i = qlen+1; i < (qlen/pn+1) * pn; ++i) _qi[i] = -1;                                   \
    }                                                                                              \
    /* for backtrack */                                                                            \
    dp_beg=abm->dp_beg, dp_end=abm->dp_end, dp_beg_sn=abm->dp_beg_sn, dp_end_sn=abm->dp_end_sn;    \
    /* index of pre-node */                                                                        \
    pre_index = (int**)_err_calloc(matrix_row_n, sizeof(int*));                                    \
    pre_n = (int*)_err_calloc(matrix_row_n, sizeof(int));                                          \
    for (index_i=beg_index+1, dp_i=1; index_i<=end_index; ++index_i, ++dp_i) {                     \
        node_id = abpoa_graph_index_to_node_id(graph, index_i);                                    \
        pre_n[dp_i] = graph->node[node_id].in_edge_n;                                              \
        pre_index[dp_i] = (int*)_err_malloc(pre_n[dp_i] * sizeof(int));                            \
        for (j = _pre_n = 0; j < pre_n[dp_i]; ++j) {                                               \
            _pre_index = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);       \
            if (index_map[_pre_index]) pre_index[dp_i][_pre_n++] = _pre_index-beg_index;           \
        }                                                                                          \
        pre_n[dp_i] = _pre_n;                                                                      \
    }                                                                                              \
}

#define simd_abpoa_free_var {                                                            \
    for (i = 0; i < matrix_row_n; ++i) free(pre_index[i]); free(pre_index); free(pre_n); \
    SIMDFree(PRE_MASK); SIMDFree(SUF_MIN); SIMDFree(PRE_MIN);                            \
}                                                                                        \

#define simd_abpoa_lg_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd) \
    simd_abpoa_var(score_t, sp, SIMDSetOne, SIMDShiftOneN);                \
    simd_abpoa_lg_only_var(score_t, SIMDSetOne, SIMDAdd);                  \
    simd_abpoa_init_var(score_t);

#define simd_abpoa_ag_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd) \
    simd_abpoa_var(score_t, sp, SIMDSetOne, SIMDShiftOneN);                \
    simd_abpoa_ag_only_var(score_t, SIMDSetOne, SIMDAdd);                  \
    simd_abpoa_init_var(score_t);

#define simd_abpoa_cg_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd) \
    simd_abpoa_var(score_t, sp, SIMDSetOne, SIMDShiftOneN);                \
    simd_abpoa_cg_only_var(score_t, SIMDSetOne, SIMDAdd);                  \
    simd_abpoa_init_var(score_t);

#define simd_abpoa_lg_first_row {                                                                       \
    /* fill the first row */                                                                            \
    if (abpt->wb >= 0) {                                                                                \
        graph->node_id_to_max_pos_left[beg_node_id] = graph->node_id_to_max_pos_right[beg_node_id] = 0; \
        for (i = 0; i < graph->node[beg_node_id].out_edge_n; ++i) { /* set max pos for out_id */        \
            int out_id = graph->node[beg_node_id].out_id[i];                                            \
            if (index_map[abpoa_graph_node_id_to_index(graph, out_id)])                                 \
                graph->node_id_to_max_pos_left[out_id] = graph->node_id_to_max_pos_right[out_id] = 1;   \
        }                                                                                               \
        dp_beg[0] = 0, dp_end[0] = GET_AD_DP_END(graph, w, beg_node_id, end_node_id, qlen);             \
    } else {                                                                                            \
        dp_beg[0] = 0, dp_end[0] = qlen;                                                                \
    }                                                                                                   \
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;                                       \
    dp_h = DP_H; _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);                                         \
}

#define simd_abpoa_ag_first_row {                                                                       \
    /* fill the first row */                                                                            \
    if (abpt->wb >= 0) {                                                                                \
        graph->node_id_to_max_pos_left[beg_node_id] = graph->node_id_to_max_pos_right[beg_node_id] = 0; \
        for (i = 0; i < graph->node[beg_node_id].out_edge_n; ++i) { /* set max pos for out_id */        \
            int out_id = graph->node[beg_node_id].out_id[i];                                            \
            if (index_map[abpoa_graph_node_id_to_index(graph, out_id)])                                 \
                graph->node_id_to_max_pos_left[out_id] = graph->node_id_to_max_pos_right[out_id] = 1;   \
        }                                                                                               \
        dp_beg[0] = 0, dp_end[0] = GET_AD_DP_END(graph, w, beg_node_id, end_node_id, qlen);             \
    } else {                                                                                            \
        dp_beg[0] = 0, dp_end[0] = qlen;                                                                \
    }                                                                                                   \
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;                                       \
    dp_h = DP_HEF; dp_e1 = dp_h + dp_sn; dp_f1 = dp_e1 + dp_sn;                                         \
    _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);                                                      \
}

#define simd_abpoa_cg_first_row {                                                                       \
    /* fill the first row */                                                                            \
    if (abpt->wb >= 0) {                                                                                \
        graph->node_id_to_max_pos_left[beg_node_id] = graph->node_id_to_max_pos_right[beg_node_id] = 0; \
        for (i = 0; i < graph->node[beg_node_id].out_edge_n; ++i) { /* set max pos for out_id */        \
            int out_id = graph->node[beg_node_id].out_id[i];                                            \
            if (index_map[abpoa_graph_node_id_to_index(graph, out_id)])                                 \
                graph->node_id_to_max_pos_left[out_id] = graph->node_id_to_max_pos_right[out_id] = 1;   \
        }                                                                                               \
        dp_beg[0] = 0, dp_end[0] = GET_AD_DP_END(graph, w, beg_node_id, end_node_id, qlen);             \
    } else {                                                                                            \
        dp_beg[0] = 0, dp_end[0] = qlen;                                                                \
    }                                                                                                   \
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;                                       \
    dp_h = DP_H2E2F; dp_e1 = dp_h+dp_sn; dp_e2 = dp_e1+dp_sn; dp_f1 = dp_e2+dp_sn; dp_f2 = dp_f1+dp_sn; \
    _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);                                                      \
}

#define simd_abpoa_lg_first_dp(score_t) {                                   \
    simd_abpoa_lg_first_row;                                                \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                             \
        for (i = 0; i <= _end_sn; ++i)                                      \
            dp_h[i] = zero;                                                 \
    } else {                                                                \
        for (i = 0; i <= _end_sn; ++i) {                                    \
            dp_h[i] = SIMD_INF_MIN;                                         \
        }                                                                   \
        _dp_h = (score_t*)dp_h;                                             \
        for (i = 0; i <= dp_end[0]; ++i) { /* no SIMD parallelization */    \
            _dp_h[i] = -gap_ext1 * i;                                       \
        }                                                                   \
    }                                                                       \
}

#define simd_abpoa_ag_first_dp(score_t) {                                   \
    simd_abpoa_ag_first_row;                                                \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                             \
        for (i = 0; i <= _end_sn; ++i)                                      \
            dp_h[i] = dp_e1[i] = dp_f1[i] = zero;                           \
    } else {                                                                \
        for (i = 0; i <= _end_sn; ++i) {                                    \
            dp_h[i] = SIMD_INF_MIN; dp_e1[i] = SIMD_INF_MIN;                \
        }                                                                   \
        _dp_h=(score_t*)dp_h,_dp_e1=(score_t*)dp_e1,_dp_f1=(score_t*)dp_f1; \
        _dp_h[0] = 0; _dp_e1[0] = -(gap_oe1), _dp_f1[0] = inf_min;          \
        for (i = 1; i <= dp_end[0]; ++i) { /* no SIMD parallelization */    \
            _dp_f1[i] = -gap_open1 - gap_ext1 * i;                          \
            _dp_h[i] = -gap_open1 - gap_ext1 * i;                           \
        }                                                                   \
    }                                                                       \
}

#define simd_abpoa_cg_first_dp(score_t) {                                             \
    simd_abpoa_cg_first_row;                                                          \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                       \
        for (i = 0; i <= _end_sn; ++i)                                                \
            dp_h[i] = dp_e1[i] = dp_e2[i] = dp_f1[i] = dp_f2[i] = zero;               \
    } else {                                                                          \
        for (i = 0; i <= _end_sn; ++i) {                                              \
            dp_h[i] = SIMD_INF_MIN; dp_e1[i] = SIMD_INF_MIN; dp_e2[i] = SIMD_INF_MIN; \
        }                                                                             \
        _dp_h = (score_t*)dp_h, _dp_e1 = (score_t*)dp_e1, _dp_e2 = (score_t*)dp_e2;   \
        _dp_f1 = (score_t*)dp_f1, _dp_f2 = (score_t*)dp_f2;                           \
        _dp_h[0] = 0; _dp_e1[0] = -(gap_oe1); _dp_e2[0] = -(gap_oe2);                 \
        _dp_f1[0] = _dp_f2[0] = inf_min;                                              \
        for (i = 1; i <= dp_end[0]; ++i) { /* no SIMD parallelization */              \
            _dp_f1[i] = -gap_open1 - gap_ext1 * i;                                    \
            _dp_f2[i] = -gap_open2 - gap_ext2 * i;                                    \
            _dp_h[i] = MAX_OF_TWO(_dp_f1[i], _dp_f2[i]);                              \
        }                                                                             \
    }                                                                                 \
}

// mask[pn], suf_min[pn], pre_min[logN]
#define SIMD_SET_F(F, log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN) { \
    if (set_num == pn) {                                                                                               \
        F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[0]), SIMDShiftOneN), PRE_MIN[1]));                     \
        if (log_n > 1) {                                                                                               \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[1]), SIMDShiftOneN<<1), PRE_MIN[2]));              \
        } if (log_n > 2) {                                                                                             \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[2]), SIMDShiftOneN<<2), PRE_MIN[4]));              \
        } if (log_n > 3) {                                                                                             \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[3]), SIMDShiftOneN<<3), PRE_MIN[8]));              \
        } if (log_n > 4) {                                                                                             \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[4]), SIMDShiftOneN<<4), PRE_MIN[16]));             \
        } if (log_n > 5) {                                                                                             \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[5]), SIMDShiftOneN<<5), PRE_MIN[32]));             \
        }                                                                                                              \
    } else { /*suffix MIN_INF*/                                                                                                                                    \
        int cov_bit = set_num;                                                                                                                                     \
        F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[0]), SIMDShiftOneN), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[1])));         \
        if (log_n > 1) {                                                                                                                                           \
            cov_bit += 2;                                                                                                                                          \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[1]), SIMDShiftOneN<<1), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[2])));  \
        } if (log_n > 2) {                                                                                                                                         \
            cov_bit += 4;                                                                                                                                          \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[2]), SIMDShiftOneN<<2), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[4])));  \
        } if (log_n > 3) {                                                                                                                                         \
            cov_bit += 8;                                                                                                                                          \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[3]), SIMDShiftOneN<<3), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[8])));  \
        } if (log_n > 4) {                                                                                                                                         \
            cov_bit += 16;                                                                                                                                         \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[4]), SIMDShiftOneN<<4), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[16]))); \
        } if (log_n > 5) {                                                                                                                                         \
            cov_bit += 32;                                                                                                                                         \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[5]), SIMDShiftOneN<<5), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[32]))); \
        }                                                                                                                                                          \
    }                                                                                                                                                              \
}

#define simd_abpoa_lg_dp(score_t, SIMDSetOne, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub) {                                         \
    node_id = abpoa_graph_index_to_node_id(graph, index_i);                                                                       \
    SIMDi *q = qp + graph->node[node_id].base * qp_sn, first, remain;                                                             \
    dp_h = &DP_H[dp_i * dp_sn]; _dp_h = (score_t*)dp_h;                                                                           \
    int min_pre_beg, min_pre_beg_sn, max_pre_end_sn; int path_score = 0;                                                          \
    if (abpt->wb < 0) {                                                                                                           \
        beg = dp_beg[dp_i] = 0, end = dp_end[dp_i] = qlen;                                                                        \
        beg_sn = dp_beg_sn[dp_i] = (dp_beg[dp_i])/pn; end_sn = dp_end_sn[dp_i] = (dp_end[dp_i])/pn;                               \
        min_pre_beg_sn = 0, max_pre_end_sn = end_sn;                                                                              \
    } else {                                                                                                                      \
        beg = GET_AD_DP_BEGIN(graph, w, node_id, end_node_id, qlen), end = GET_AD_DP_END(graph, w, node_id, end_node_id, qlen);   \
        beg_sn = beg / pn; min_pre_beg = min_pre_beg_sn = INT32_MAX, max_pre_end_sn = -1;                                         \
        for (i = 0; i < pre_n[dp_i]; ++i) {                                                                                       \
            pre_i = pre_index[dp_i][i];                                                                                           \
            if (min_pre_beg > dp_beg[pre_i]) min_pre_beg = dp_beg[pre_i], min_pre_beg_sn = dp_beg_sn[pre_i];                      \
            if (max_pre_end_sn < dp_end_sn[pre_i]) max_pre_end_sn = dp_end_sn[pre_i];                                             \
        }                                                                                                                         \
        if (beg_sn < min_pre_beg_sn) {                                                                                            \
            beg = min_pre_beg; beg_sn = min_pre_beg_sn;                                                                           \
        }                                                                                                                         \
        dp_beg_sn[dp_i] = beg_sn; dp_beg[dp_i] = beg; end_sn = dp_end_sn[dp_i] = end/pn; dp_end[dp_i] = end;                      \
    }                                                                                                                             \
    /* loop query */                                                                                                              \
    /* first pre_node */                                                                                                          \
    pre_i = pre_index[dp_i][0];                                                                                                   \
    if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, node_id, 0);                                         \
    SIMDi dp_path_score = SIMDSetOne(path_score);                                                                                 \
    pre_dp_h = DP_H + pre_i * dp_sn;                                                                                              \
    pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i];                                                                       \
    /* set M from (pre_i, q_i-1), E from (pre_i, q_i) */                                                                          \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                                   \
        _beg_sn = 0, _end_sn = end_sn; first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                \
    } else {                                                                                                                      \
        if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);      \
        else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                            \
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                                  \
        for (i = beg_sn; i < _beg_sn; ++i) dp_h[i] = SIMD_INF_MIN;                                                                \
        for (i = _end_sn+1; i <= MIN_OF_TWO(end_sn+1, dp_sn-1); ++i) dp_h[i] = SIMD_INF_MIN;                                      \
    }                                                                                                                             \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                    \
        remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                    \
        dp_h[sn_i] = SIMDMax(SIMDAdd(SIMDOri(first, remain), SIMDAdd(dp_path_score, q[sn_i])), SIMDAdd(dp_path_score, SIMDSub(pre_dp_h[sn_i], GAP_E1))); \
        first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                     \
    }                                                                                                                             \
    /* get max m and e */                                                                                                         \
    for (i = 1; i < pre_n[dp_i]; ++i) {                                                                                           \
        pre_i = pre_index[dp_i][i];                                                                                               \
        if (abpt->inc_path_score) {                                                                                               \
            path_score = abpoa_get_incre_path_score(graph, node_id, i);                                                           \
            dp_path_score = SIMDSetOne(path_score);                                                                               \
        }                                                                                                                         \
        pre_dp_h = DP_H + pre_i * dp_sn;                                                                                          \
        pre_end = dp_end[pre_i];                                                                                                  \
        pre_beg_sn = dp_beg_sn[pre_i];                                                                                            \
        /* set M from (pre_i, q_i-1), E from (pre_i, q_i) */                                                                      \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                                           \
        } else {                                                                                                                  \
            if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);  \
            else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                        \
            _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                              \
        }                                                                                                                         \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                \
            remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                \
            dp_h[sn_i] = SIMDMax(SIMDAdd(SIMDOri(first, remain), SIMDAdd(dp_path_score, q[sn_i])), SIMDMax(SIMDAdd(dp_path_score, SIMDSub(pre_dp_h[sn_i], GAP_E1)), dp_h[sn_i])); \
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                 \
        } /* now we have max(h,e) stored at dp_h */                                                                               \
    }                                                                                                                             \
    /* new F start */                                                                                                             \
    for (i = beg_sn * pn; i < beg; ++i) _dp_h[i] = inf_min;                                                                       \
    for (i = end+1; i < (end_sn+1)*pn; ++i) _dp_h[i] = inf_min;                                                                   \
    first = SIMDOri(SIMDAndi(dp_h[beg_sn], PRE_MASK[0]), SUF_MIN[0]);                                                             \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) {                                                                                 \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            set_num = pn;                                                                                                         \
        } else {                                                                                                                  \
            if (sn_i < min_pre_beg_sn) {                                                                                          \
                _err_fatal_simple(__func__, "sn_i < min_pre_beg_sn\n");                                                           \
            } else if (sn_i > max_pre_end_sn) {                                                                                   \
                set_num = sn_i == max_pre_end_sn+1 ? 1 : 0;                                                                       \
            } else set_num = pn;                                                                                                  \
        }                                                                                                                         \
        dp_h[sn_i] = SIMDMax(dp_h[sn_i], first);                                                                                  \
        if (sn_i == end_sn) for (i = end+1; i < (end_sn+1)*pn; ++i) _dp_h[i] = inf_min;                                           \
        SIMD_SET_F(dp_h[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN);    \
        first = SIMDOri(SIMDAndi(SIMDShiftRight(SIMDSub(dp_h[sn_i], GAP_E1), SIMDTotalBytes-SIMDShiftOneN), PRE_MASK[0]), SUF_MIN[0]); \
    }                                                                                                                             \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) for (sn_i = 0; sn_i <= end_sn; ++sn_i) dp_h[sn_i] = SIMDMax(zero, dp_h[sn_i]);      \
}

#define simd_abpoa_ag_dp(score_t, SIMDSetOne, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater, SIMDSetIfEqual) { \
    node_id = abpoa_graph_index_to_node_id(graph, index_i);                                                                       \
    SIMDi *q = qp + graph->node[node_id].base * qp_sn, first, remain;                                                             \
    dp_h = DP_HEF + dp_i * 3 * dp_sn; dp_e1 = dp_h + dp_sn; dp_f1 = dp_e1 + dp_sn;                                                \
    _dp_h = (score_t*)dp_h, _dp_e1 = (score_t*)dp_e1, _dp_f1 = (score_t*)dp_f1;                                                   \
    int min_pre_beg, min_pre_beg_sn, max_pre_end_sn, path_score = 0;                                                              \
    if (abpt->wb < 0) {                                                                                                           \
        beg = dp_beg[dp_i] = 0, end = dp_end[dp_i] = qlen;                                                                        \
        beg_sn = dp_beg_sn[dp_i] = (dp_beg[dp_i])/pn; end_sn = dp_end_sn[dp_i] = (dp_end[dp_i])/pn;                               \
        min_pre_beg_sn = 0, max_pre_end_sn = end_sn;                                                                              \
    } else {                                                                                                                      \
        beg = GET_AD_DP_BEGIN(graph, w, node_id, end_node_id, qlen), end = GET_AD_DP_END(graph, w, node_id, end_node_id, qlen);   \
        beg_sn = beg / pn; min_pre_beg = min_pre_beg_sn = INT32_MAX, max_pre_end_sn = -1;                                         \
        for (i = 0; i < pre_n[dp_i]; ++i) {                                                                                       \
            pre_i = pre_index[dp_i][i];                                                                                           \
            if (min_pre_beg > dp_beg[pre_i]) min_pre_beg = dp_beg[pre_i], min_pre_beg_sn = dp_beg_sn[pre_i];                      \
            if (max_pre_end_sn < dp_end_sn[pre_i]) max_pre_end_sn = dp_end_sn[pre_i];                                             \
        }                                                                                                                         \
        if (beg_sn < min_pre_beg_sn) {                                                                                            \
            beg_sn = min_pre_beg_sn; beg = min_pre_beg;                                                                           \
        }                                                                                                                         \
        dp_beg_sn[dp_i] = beg_sn; dp_beg[dp_i] = beg; end_sn = dp_end_sn[dp_i] = end/pn; dp_end[dp_i] = end;                      \
    }                                                                                                                             \
    /* loop query */                                                                                                              \
    /* first pre_node */                                                                                                          \
    pre_i = pre_index[dp_i][0];                                                                                                   \
    if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, node_id, 0);                                         \
    SIMDi dp_path_score = SIMDSetOne(path_score);                                                                                 \
    pre_dp_h = DP_HEF + pre_i * 3 * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn;                                                          \
    pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                        \
    /* set M from (pre_i, q_i-1) */                                                                                               \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                                   \
        _beg_sn = 0, _end_sn = end_sn; first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                \
    } else {                                                                                                                      \
        if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);      \
        else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                            \
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                                  \
        for (i = beg_sn; i < _beg_sn; ++i) dp_h[i] = SIMD_INF_MIN;                                                                \
        for (i = _end_sn+1; i <= MIN_OF_TWO(end_sn+1, dp_sn-1); ++i) dp_h[i] = SIMD_INF_MIN;                                      \
    }                                                                                                                             \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                    \
        remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                    \
        dp_h[sn_i] = SIMDAdd(SIMDOri(first, remain), dp_path_score);                                                              \
        first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                     \
    }                                                                                                                             \
    /* set E from (pre_i, q_i) */                                                                                                 \
    if (abpt->align_mode != ABPOA_LOCAL_MODE) {                                                                                   \
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);                                                                                 \
        for (i = beg_sn; i < _beg_sn; ++i) dp_e1[i] = SIMD_INF_MIN;                                                               \
        for (i = _end_sn+1; i <= end_sn; ++i) dp_e1[i] = SIMD_INF_MIN;                                                            \
    }                                                                                                                             \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i)   /* SIMD parallelization */                                                    \
        dp_e1[sn_i] = SIMDAdd(pre_dp_e1[sn_i], dp_path_score);                                                                    \
    /* get max m and e */                                                                                                         \
    for (i = 1; i < pre_n[dp_i]; ++i) {                                                                                           \
        pre_i = pre_index[dp_i][i];                                                                                               \
        if (abpt->inc_path_score) {                                                                                               \
            path_score = abpoa_get_incre_path_score(graph, node_id, i);                                                           \
            dp_path_score = SIMDSetOne(path_score);                                                                               \
        }                                                                                                                         \
        pre_dp_h = DP_HEF + pre_i * 3 * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn;                                                      \
        pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                    \
        /* set M from (pre_i, q_i-1) */                                                                                           \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                                           \
        } else {                                                                                                                  \
            if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);  \
            else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                        \
            _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                              \
        }                                                                                                                         \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                \
            remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                \
            dp_h[sn_i] = SIMDMax(SIMDAdd(SIMDOri(first, remain), dp_path_score), dp_h[sn_i]);                                     \
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                 \
        }                                                                                                                         \
        /* set E from (pre_i, q_i) */                                                                                             \
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);                                                                                 \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i)   /* SIMD parallelization */                                                \
            dp_e1[sn_i] = SIMDMax(SIMDAdd(pre_dp_e1[sn_i], dp_path_score), dp_e1[sn_i]);                                          \
    }                                                                                                                             \
    /* compare M, E, and F */                                                                                                     \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */                                                      \
        dp_h[sn_i] = SIMDAdd(dp_h[sn_i], q[sn_i]);                                                                                \
    }                                                                                                                             \
    for (i = beg_sn * pn; i < beg; ++i) _dp_h[i] = _dp_e1[i] = inf_min;                                                           \
    for (i = end+1; i < (end_sn+1)*pn; ++i) _dp_h[i] = _dp_e1[i] = inf_min;                                                       \
    /* new F start */                                                                                                             \
    first = SIMDShiftRight(SIMDShiftLeft(dp_h[beg_sn], SIMDTotalBytes-SIMDShiftOneN), SIMDTotalBytes-SIMDShiftOneN);              \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) {                                                                                 \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            set_num  = pn;                                                                                                        \
        } else {                                                                                                                  \
            if (sn_i < min_pre_beg_sn) {                                                                                          \
                _err_fatal_simple(__func__, "sn_i < min_pre_beg_sn\n");                                                           \
            } else if (sn_i > max_pre_end_sn) {                                                                                   \
                set_num = sn_i == max_pre_end_sn+1 ? 2 : 1;                                                                       \
            } else set_num = pn;                                                                                                  \
        }                                                                                                                         \
        /* F = (H << 1 | x) - OE */                                                                                               \
        dp_f1[sn_i] = SIMDSub(SIMDOri(SIMDShiftLeft(dp_h[sn_i], SIMDShiftOneN), first), GAP_OE1);                                 \
        /* F = max{F, (F-e)<<1}, F = max{F, (F-2e)<<2} ... */                                                                     \
        SIMD_SET_F(dp_f1[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN);   \
        /* x = max{H, F+o} */                                                                                                     \
        first = SIMDShiftRight(SIMDMax(dp_h[sn_i], SIMDAdd(dp_f1[sn_i], GAP_O1)), SIMDTotalBytes-SIMDShiftOneN);                  \
        /* H = max{H, F} */                                                                                                       \
        dp_h[sn_i] = SIMDMax(dp_h[sn_i], dp_e1[sn_i]); SIMDi tmp = dp_h[sn_i];                                                    \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            dp_h[sn_i] = SIMDMax(zero, SIMDMax(dp_h[sn_i], dp_f1[sn_i]));                                                         \
            if (sn_i == end_sn) for (i = end+1; i < (end_sn+1)*pn; ++i) _dp_h[i] = _dp_e1[i] = inf_min;                           \
            SIMDSetIfEqual(dp_e1[sn_i], dp_h[sn_i],tmp, SIMDMax(SIMDSub(dp_e1[sn_i],GAP_E1), SIMDSub(dp_h[sn_i],GAP_OE1)),zero);  \
        } else {                                                                                                                  \
            dp_h[sn_i] = SIMDMax(dp_h[sn_i], dp_f1[sn_i]);                                                                        \
            if (sn_i == end_sn) for (i = end+1; i < (end_sn+1)*pn; ++i) _dp_h[i] = _dp_e1[i] = inf_min;                           \
            SIMDSetIfEqual(dp_e1[sn_i], dp_h[sn_i],tmp, SIMDMax(SIMDSub(dp_e1[sn_i],GAP_E1), SIMDSub(dp_h[sn_i],GAP_OE1)),SIMD_INF_MIN); \
        }                                                                                                                         \
    }                                                                                                                             \
}

#define simd_abpoa_cg_dp(score_t, SIMDSetOne, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater, SIMDSetIfEqual) { \
    node_id = abpoa_graph_index_to_node_id(graph, index_i);                                                                       \
    SIMDi *q = qp + graph->node[node_id].base * qp_sn, first, remain;                                                             \
    dp_h = DP_H2E2F+dp_i*5*dp_sn; dp_e1 = dp_h+dp_sn; dp_e2 = dp_e1+dp_sn; dp_f1 = dp_e2+dp_sn; dp_f2 = dp_f1+dp_sn;              \
    _dp_h=(score_t*)dp_h, _dp_e1=(score_t*)dp_e1, _dp_e2=(score_t*)dp_e2, _dp_f1=(score_t*)dp_f1, _dp_f2=(score_t*)dp_f2;         \
    int min_pre_beg, min_pre_beg_sn, max_pre_end_sn, path_score = 0;                                                              \
    if (abpt->wb < 0) {                                                                                                           \
        beg = dp_beg[dp_i] = 0, end = dp_end[dp_i] = qlen;                                                                        \
        beg_sn = dp_beg_sn[dp_i] = beg/pn; end_sn = dp_end_sn[dp_i] = end/pn;                                                     \
        min_pre_beg_sn = 0, max_pre_end_sn = end_sn;                                                                              \
    } else {                                                                                                                      \
        beg = GET_AD_DP_BEGIN(graph, w, node_id, end_node_id, qlen), end = GET_AD_DP_END(graph, w, node_id, end_node_id, qlen);   \
        if (abpt->verbose >= ABPOA_LONG_DEBUG_VERBOSE)                                                                            \
            fprintf(stderr, "index: %d (node: %d): raw beg: %d, end: %d\n", index_i-beg_index, node_id, beg, end);                \
        beg_sn = beg / pn; min_pre_beg = min_pre_beg_sn = INT32_MAX, max_pre_end_sn = -1;                                         \
        for (i = 0; i < pre_n[dp_i]; ++i) {                                                                                       \
            pre_i = pre_index[dp_i][i];                                                                                           \
            if (abpt->verbose >= ABPOA_LONG_DEBUG_VERBOSE)                                                                        \
                fprintf(stderr, "\tpre_i: %d, pre_dp_beg: %d, pre_dp_beg_sn: %d\n", pre_i, dp_beg[pre_i], dp_beg_sn[pre_i]);      \
            if (min_pre_beg > dp_beg[pre_i]) min_pre_beg = dp_beg[pre_i], min_pre_beg_sn = dp_beg_sn[pre_i];                      \
            if (max_pre_end_sn < dp_end_sn[pre_i]) max_pre_end_sn = dp_end_sn[pre_i];                                             \
        }                                                                                                                         \
        if (beg_sn < min_pre_beg_sn) {                                                                                            \
            assert(beg < min_pre_beg); beg = min_pre_beg; beg_sn = min_pre_beg_sn;                                                \
        }                                                                                                                         \
        dp_beg_sn[dp_i] = beg_sn; dp_beg[dp_i] = beg; end_sn = dp_end_sn[dp_i] = end/pn; dp_end[dp_i] = end;                      \
        if (abpt->verbose >= ABPOA_LONG_DEBUG_VERBOSE) {                                                                          \
            fprintf(stderr, "index: %d (node: %d): beg: %d, end: %d, beg_sn: %d, end_sn: %d, max_left: %d, max_right: %d\n", index_i-beg_index, node_id, beg, end, beg_sn, end_sn, abpoa_graph_node_id_to_max_pos_left(graph, node_id), abpoa_graph_node_id_to_max_pos_right(graph, node_id)); \
        }                                                                                                                         \
    }                                                                                                                             \
    /* loop query */                                                                                                              \
    /* first pre_node */                                                                                                          \
    pre_i = pre_index[dp_i][0];                                                                                                   \
    if (abpt->inc_path_score) path_score = abpoa_get_incre_path_score(graph, node_id, 0);                                         \
    SIMDi dp_path_score = SIMDSetOne(path_score);                                                                                 \
    pre_dp_h = DP_H2E2F + pre_i * 5 * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn; pre_dp_e2 = pre_dp_e1 + dp_sn;                         \
    pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                        \
    /* set M from (pre_i, q_i-1) */                                                                                               \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                                   \
        _beg_sn = 0, _end_sn = end_sn; first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                \
    } else {                                                                                                                      \
        if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);      \
        else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                            \
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                                  \
        for (i = beg_sn; i < _beg_sn; ++i) dp_h[i] = SIMD_INF_MIN;                                                                \
        for (i = _end_sn+1; i <= MIN_OF_TWO(end_sn+1, dp_sn-1); ++i) dp_h[i] = SIMD_INF_MIN;                                      \
    }                                                                                                                             \
 /* fprintf(stderr, "1 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                     \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                    \
        remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                    \
        dp_h[sn_i] = SIMDAdd(SIMDOri(first, remain), dp_path_score);                                                              \
        first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                     \
    }                                                                                                                             \
    /* set E from (pre_i, q_i) */                                                                                                 \
    if (abpt->align_mode != ABPOA_LOCAL_MODE) {                                                                                   \
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);                                                                                 \
        for (i = beg_sn; i < _beg_sn; ++i) dp_e1[i] = SIMD_INF_MIN, dp_e2[i] = SIMD_INF_MIN;                                      \
        for (i = _end_sn+1; i <= end_sn; ++i) dp_e1[i] = SIMD_INF_MIN, dp_e2[i] = SIMD_INF_MIN;                                   \
    }                                                                                                                             \
 /* fprintf(stderr, "2 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                     \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                    \
        dp_e1[sn_i] = SIMDAdd(pre_dp_e1[sn_i], dp_path_score);                                                                    \
        dp_e2[sn_i] = SIMDAdd(pre_dp_e2[sn_i], dp_path_score);                                                                    \
    }                                                                                                                             \
    /* get max m and e */                                                                                                         \
    for (i = 1; i < pre_n[dp_i]; ++i) {                                                                                           \
        pre_i = pre_index[dp_i][i];                                                                                               \
        if (abpt->inc_path_score) {                                                                                               \
            path_score = abpoa_get_incre_path_score(graph, node_id, i);                                                           \
            dp_path_score = SIMDSetOne(path_score);                                                                               \
        }                                                                                                                         \
        pre_dp_h = DP_H2E2F + (pre_i * 5) * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn; pre_dp_e2 = pre_dp_e1 + dp_sn;                   \
        pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                    \
        /* set M from (pre_i, q_i-1) */                                                                                           \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                                           \
        } else {                                                                                                                  \
            if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);  \
            else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                        \
            _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                              \
        }                                                                                                                         \
     /* fprintf(stderr, "3 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                 \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                \
            remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                \
            dp_h[sn_i] = SIMDMax(SIMDAdd(SIMDOri(first, remain), dp_path_score), dp_h[sn_i]);                                     \
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                 \
        }                                                                                                                         \
        /* set E from (pre_i, q_i) */                                                                                             \
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);                                                                                 \
     /* fprintf(stderr, "4 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                 \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                \
            dp_e1[sn_i] = SIMDMax(SIMDAdd(pre_dp_e1[sn_i], dp_path_score), dp_e1[sn_i]);                                          \
            dp_e2[sn_i] = SIMDMax(SIMDAdd(pre_dp_e2[sn_i], dp_path_score), dp_e2[sn_i]);                                          \
        }                                                                                                                         \
    }                                                                                                                             \
    /* compare M, E, and F */                                                                                                     \
 /* fprintf(stderr, "5 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, beg_sn, end_sn); */                                       \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */                                                      \
        dp_h[sn_i] = SIMDAdd(dp_h[sn_i], q[sn_i]);                                                                                \
    }                                                                                                                             \
    for (i = beg_sn * pn; i < beg; ++i) _dp_h[i] = _dp_e1[i] = _dp_e2[i] = inf_min;                                               \
    for (i = end+1; i < (end_sn+1)*pn; ++i) _dp_h[i] = _dp_e1[i] = _dp_e2[i] = inf_min;                                           \
    /* new F start */                                                                                                             \
    first = SIMDShiftRight(SIMDShiftLeft(dp_h[beg_sn], SIMDTotalBytes-SIMDShiftOneN), SIMDTotalBytes-SIMDShiftOneN);              \
    SIMDi first2 = first;                                                                                                         \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) {                                                                                 \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) set_num = pn;                                                                   \
        else {                                                                                                                    \
            if (sn_i < min_pre_beg_sn) {                                                                                          \
                _err_fatal_simple(__func__, "sn_i < min_pre_beg_sn\n");                                                           \
            } else if (sn_i > max_pre_end_sn) {                                                                                   \
                set_num = sn_i == max_pre_end_sn+1 ? 2 : 1;                                                                       \
            } else set_num = pn;                                                                                                  \
        }                                                                                                                         \
        /* H = max{H, E} */                                                                                                       \
        dp_h[sn_i] =  SIMDMax(SIMDMax(dp_h[sn_i], dp_e1[sn_i]), dp_e2[sn_i]);                                                     \
        /* F = (H << 1 | x) - OE */                                                                                               \
        dp_f1[sn_i] = SIMDSub(SIMDOri(SIMDShiftLeft(dp_h[sn_i], SIMDShiftOneN), first), GAP_OE1);                                 \
        dp_f2[sn_i] = SIMDSub(SIMDOri(SIMDShiftLeft(dp_h[sn_i], SIMDShiftOneN), first2), GAP_OE2);                                \
        /* F = max{F, (F-e)<<1}, F = max{F, (F-2e)<<2} ... */                                                                     \
        SIMD_SET_F(dp_f1[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN);   \
        SIMD_SET_F(dp_f2[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E2S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN);   \
        /* x = max{H, F+o} */                                                                                                     \
        first = SIMDShiftRight(SIMDMax(dp_h[sn_i], SIMDAdd(dp_f1[sn_i], GAP_O1)), SIMDTotalBytes-SIMDShiftOneN);                  \
        first2 = SIMDShiftRight(SIMDMax(dp_h[sn_i], SIMDAdd(dp_f2[sn_i], GAP_O2)), SIMDTotalBytes-SIMDShiftOneN);                 \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            dp_h[sn_i] = SIMDMax(zero, SIMDMax(dp_h[sn_i], SIMDMax(dp_f1[sn_i], dp_f2[sn_i])));                                   \
            if (sn_i == end_sn) for (i = end+1; i < (end_sn+1)*pn; ++i) _dp_h[i] = _dp_e1[i] = _dp_e2[i] = inf_min;               \
            dp_e1[sn_i] = SIMDMax(zero,SIMDMax(SIMDSub(dp_e1[sn_i],GAP_E1),SIMDSub(dp_h[sn_i],GAP_OE1)));                         \
            dp_e2[sn_i] = SIMDMax(zero,SIMDMax(SIMDSub(dp_e2[sn_i],GAP_E2),SIMDSub(dp_h[sn_i],GAP_OE2)));                         \
        } else {                                                                                                                  \
            /* H = max{H, F}    */                                                                                                \
            dp_h[sn_i] = SIMDMax(dp_h[sn_i], SIMDMax(dp_f1[sn_i], dp_f2[sn_i]));                                                  \
            if (sn_i == end_sn) for (i = end+1; i < (end_sn+1)*pn; ++i) _dp_h[i] = _dp_e1[i] = _dp_e2[i] = inf_min;               \
            /* e for next cell */                                                                                                 \
            dp_e1[sn_i] = SIMDMax(SIMDSub(dp_e1[sn_i],GAP_E1),SIMDSub(dp_h[sn_i],GAP_OE1));                                       \
            dp_e2[sn_i] = SIMDMax(SIMDSub(dp_e2[sn_i],GAP_E2),SIMDSub(dp_h[sn_i],GAP_OE2));                                       \
        }                                                                                                                         \
    }                                                                                                                             \
}

#define set_global_max_score(score, i, j) {         \
    if (score > best_score) {                       \
        best_score = score; best_i = i; best_j = j; \
    }                                               \
}

#define set_extend_max_score(score, i, j) {                                                              \
    if (score > best_score) {                                                                            \
        best_score = score; best_i = i; best_j = j; best_id = node_id;                                   \
    } else if (abpt->zdrop > 0) {                                                                        \
        int delta_index = graph->node_id_to_max_remain[best_id] - graph->node_id_to_max_remain[node_id]; \
        if (best_score - score > abpt->zdrop + gap_ext1 * abs(delta_index-(j-best_j)))                   \
            break;                                                                                       \
    }                                                                                                    \
}

#define simd_abpoa_global_get_max(score_t, DP_M, dp_sn) {      \
    int end, in_id, in_index, in_dp_i;                         \
    for (i = 0; i < graph->node[end_node_id].in_edge_n; ++i) { \
        in_id = graph->node[end_node_id].in_id[i];             \
        in_index = abpoa_graph_node_id_to_index(graph, in_id); \
        if (index_map[in_index] == 0) continue;                \
        in_dp_i = in_index - beg_index;                        \
        dp_h = DP_M + in_dp_i * dp_sn;                         \
        _dp_h = (score_t*)dp_h;                                \
        if (qlen > dp_end[in_dp_i]) end = dp_end[in_dp_i];     \
        else end = qlen;                                       \
        set_global_max_score(_dp_h[end], in_dp_i, end);        \
    }                                                          \
}

#define simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater) { \
    /* select max dp_h */                                                    \
    max = inf_min, left_max_i = right_max_i = -1;                            \
    _dp_h = (score_t*)dp_h, _qi = (score_t*)qi;                              \
    for (i = beg; i <= end; ++i) {                                           \
        if (_dp_h[i] > max) {                                                \
            max = _dp_h[i];                                                  \
            left_max_i = right_max_i = _qi[i];                               \
        } else if (_dp_h[i] == max) {                                        \
            right_max_i = _qi[i];                                            \
        }                                                                    \
    }                                                                        \
}

#define simd_abpoa_ada_max_i   {                                             \
    /* set max_pos_left/right for next nodes */                              \
    for (i = 0; i < graph->node[node_id].out_edge_n; ++i) {                  \
        int out_node_id = graph->node[node_id].out_id[i];                    \
        if (right_max_i+1 > graph->node_id_to_max_pos_right[out_node_id])    \
            graph->node_id_to_max_pos_right[out_node_id] = right_max_i+1;    \
        if (left_max_i+1 < graph->node_id_to_max_pos_left[out_node_id])      \
            graph->node_id_to_max_pos_left[out_node_id] = left_max_i+1;      \
    }                                                                        \
}

// TODO end_bonus for extension
// linear gap penalty: gap_open1 == 0
#define simd_abpoa_lg_align_sequence_to_graph_core(score_t, sp, SIMDSetOne, SIMDMax, SIMDAdd,   \
        SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater) {                           \
    simd_abpoa_lg_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd);                         \
    simd_abpoa_lg_first_dp(score_t);                                                            \
    for (index_i=beg_index+1, dp_i=1; index_i<end_index; ++index_i, ++dp_i) {                   \
        if (index_map[index_i] == 0) continue;                                                  \
        simd_abpoa_lg_dp(score_t, SIMDSetOne, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub);        \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                             \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                 \
            set_global_max_score(max, dp_i, left_max_i);                                        \
        }                                                                                       \
        if (abpt->align_mode == ABPOA_EXTEND_MODE) {                                            \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                 \
            set_extend_max_score(max, dp_i, right_max_i);                                       \
        }                                                                                       \
        if (abpt->wb >= 0) {                                                                    \
            if (abpt->align_mode == ABPOA_GLOBAL_MODE) {                                        \
                simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);             \
            }                                                                                   \
            simd_abpoa_ada_max_i;                                                               \
        }                                                                                       \
    }                                                                                           \
    if (abpt->align_mode == ABPOA_GLOBAL_MODE) simd_abpoa_global_get_max(score_t, DP_H, dp_sn); \
    res->best_score = best_score;                                                               \
    if (abpt->verbose >= ABPOA_DEBUG_VERBOSE) {                                                 \
        if (abpt->verbose >= ABPOA_LONG_DEBUG_VERBOSE)                                          \
            simd_abpoa_print_lg_matrix(score_t, beg_index, end_index);                          \
        fprintf(stderr, "best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);            \
    }                                                                                           \
    if (abpt->ret_cigar) simd_abpoa_lg_backtrack(score_t);                                      \
    simd_abpoa_free_var; SIMDFree(GAP_E1S);                                                     \
}

// affine gap penalty: gap_open1 > 0
#define simd_abpoa_ag_align_sequence_to_graph_core(score_t, sp, SIMDSetOne, SIMDMax, SIMDAdd,       \
        SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, SIMDSetIfEqual) {               \
    simd_abpoa_ag_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd);                             \
    simd_abpoa_ag_first_dp(score_t);                                                                \
    for (index_i=beg_index+1, dp_i=1; index_i<end_index; ++index_i, ++dp_i) {                       \
        if (index_map[index_i] == 0) continue;                                                      \
        simd_abpoa_ag_dp(score_t, SIMDSetOne, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater, SIMDSetIfEqual); \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                 \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                     \
            set_global_max_score(max, dp_i, left_max_i);                                            \
        } else if (abpt->align_mode == ABPOA_EXTEND_MODE) {                                         \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                     \
            set_extend_max_score(max, dp_i, right_max_i);                                           \
        }                                                                                           \
        if (abpt->wb >= 0) {                                                                        \
            if (abpt->align_mode == ABPOA_GLOBAL_MODE) {                                            \
                simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                 \
            }                                                                                       \
            simd_abpoa_ada_max_i;                                                                   \
        }                                                                                           \
    }                                                                                               \
    if (abpt->align_mode == ABPOA_GLOBAL_MODE) simd_abpoa_global_get_max(score_t, DP_HEF, 3*dp_sn); \
    res->best_score = best_score;                                                                   \
    if (abpt->verbose >= ABPOA_DEBUG_VERBOSE) {                                                     \
        if (abpt->verbose >= ABPOA_LONG_DEBUG_VERBOSE)                                              \
            simd_abpoa_print_ag_matrix(score_t, beg_index, end_index);                              \
        fprintf(stderr, "best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);                \
    }                                                                                               \
    if (abpt->ret_cigar) simd_abpoa_ag_backtrack(score_t);                                          \
    simd_abpoa_free_var; SIMDFree(GAP_E1S);                                                         \
}

// convex gap penalty: gap_open1 > 0 && gap_open2 > 0
#define simd_abpoa_cg_align_sequence_to_graph_core(score_t, sp, SIMDSetOne, SIMDMax, SIMDAdd,         \
        SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, SIMDSetIfEqual) {                 \
    simd_abpoa_cg_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd);                               \
    simd_abpoa_cg_first_dp(score_t);                                                                  \
    for (index_i=beg_index+1, dp_i=1; index_i<end_index; ++index_i, ++dp_i) {                         \
        if (index_map[index_i] == 0) continue;                                                        \
        simd_abpoa_cg_dp(score_t, SIMDSetOne, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater, SIMDSetIfEqual); \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                   \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                       \
            set_global_max_score(max, dp_i, left_max_i);                                              \
        } else if (abpt->align_mode == ABPOA_EXTEND_MODE) {                                           \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                       \
            set_extend_max_score(max, dp_i, right_max_i);                                             \
        }                                                                                             \
        if (abpt->wb >= 0) {                                                                          \
            if (abpt->align_mode == ABPOA_GLOBAL_MODE) {                                              \
                simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                   \
            }                                                                                         \
            simd_abpoa_ada_max_i;                                                                     \
        }                                                                                             \
    }                                                                                                 \
    if (abpt->align_mode == ABPOA_GLOBAL_MODE) simd_abpoa_global_get_max(score_t, DP_H2E2F, 5*dp_sn); \
    res->best_score = best_score;                                                                     \
    if (abpt->verbose >= ABPOA_DEBUG_VERBOSE) {                                                       \
        if (abpt->verbose >= ABPOA_LONG_DEBUG_VERBOSE)                                                \
            simd_abpoa_print_cg_matrix(score_t, beg_index, end_index);                                \
        fprintf(stderr,"best_score: (%d, %d) -> %d\n",best_i, best_j, best_score);                    \
    }                                                                                                 \
    if (abpt->ret_cigar) simd_abpoa_cg_backtrack(score_t);                                            \
    simd_abpoa_free_var; SIMDFree(GAP_E1S); SIMDFree(GAP_E2S);                                        \
}

// align query to subgraph between beg_node_id and end_node_id (both are excluded)
// generally: beg/end are the SRC/SINK_node
#ifdef ABPOA_SIMD_DISPATCH
#ifdef __AVX512BW__
#warning "Compiling with AVX512BW"
int simd_avx512_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res)
#elif defined(__AVX2__)
#warning "Compiling with AVX2"
int simd_avx2_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res)
#elif defined(__SSE4_1__)
#warning "Compiling with SSE4.1"
int simd_sse41_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res)
#elif defined(__SSE2__)
#warning "Compiling with SSE2"
int simd_sse2_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res)
#endif
#else
int simd_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res)
#endif // ~ ABPOA_SIMD_DISPATCH
{
    int old_verbose = abpt->verbose;
    if (abpt->verbose >= ABPOA_DEBUG_VERBOSE) fprintf(stderr, "beg_id: %d, end_id: %d, qlen: %d\n", beg_node_id, end_node_id, qlen);
    // if (qlen == 8957) // for debug
        // abpt->verbose = ABPOA_LONG_DEBUG_VERBOSE;
    int i, j, beg_index = ab->abg->node_id_to_index[beg_node_id], end_index = ab->abg->node_id_to_index[end_node_id];
    int gn = end_index - beg_index + 1;
    uint8_t *index_map = (uint8_t*)_err_calloc(ab->abg->node_n, sizeof(uint8_t));
    index_map[beg_index] = index_map[end_index] = 1;
    for (i = beg_index; i < end_index-1; ++i) {
        if (index_map[i] == 0) continue;
        int node_id = abpoa_graph_index_to_node_id(ab->abg, i);
        int out_n = ab->abg->node[node_id].out_edge_n;
        for (j = 0; j < out_n; ++j) {
            int out_id = ab->abg->node[node_id].out_id[j];
            index_map[abpoa_graph_node_id_to_index(ab->abg, out_id)] = 1;
        }
    }

#if __AVX512BW__
    SIMD_para_t simd_p16 = {512, 16, 5, 32, 64, -1};
    SIMD_para_t simd_p32 = {512, 32, 4, 16, 64, -1};
#elif defined(__AVX2__)
    SIMD_para_t simd_p16 = {256, 16, 4, 16, 32, -1};
    SIMD_para_t simd_p32 = {256, 32, 3,  8, 32, -1};
#else
    SIMD_para_t simd_p16 = {128, 16, 3,  8, 16, -1};
    SIMD_para_t simd_p32 = {128, 32, 2,  4, 16, -1};
#endif

    // fprintf(stderr, "reg_n: %d\n", simd_p16.reg_n);

    int32_t max_score, bits, mem_ret=0, gap_ext1 = abpt->gap_ext1, gap_ext2 = abpt->gap_ext2;
    int32_t gap_oe1 = abpt->gap_open1+gap_ext1, gap_oe2 = abpt->gap_open2+gap_ext2;
    int len = qlen > gn ? qlen : gn;
#ifdef __SIMD_DEBUG__
    simd_p32.inf_min = MAX_OF_THREE(INT32_MIN + abpt->min_mis, INT32_MIN + gap_oe1, INT32_MIN + gap_oe2) + 512 * MAX_OF_TWO(gap_ext1, gap_ext2);
    if (simd_abpoa_realloc(ab, gn, qlen, abpt, simd_p32)) return 0;
    if (abpt->gap_mode == ABPOA_CONVEX_GAP)
        abpoa_cg_global_align_sequence_to_graph_core(ab, beg_node_id, beg_index, end_node_id, end_index, index_map, qlen, query, abpt, _simd_p32, res);
#else
    max_score = MAX_OF_TWO(qlen * abpt->max_mat, len * abpt->gap_ext1 + abpt->gap_open1);
    if (max_score <= INT16_MAX - abpt->min_mis - gap_oe1 - gap_oe2) {
        simd_p16.inf_min = MAX_OF_THREE(INT16_MIN + abpt->min_mis, INT16_MIN + gap_oe1, INT16_MIN + gap_oe2) + 512 * MAX_OF_TWO(gap_ext1, gap_ext2);
        mem_ret = simd_abpoa_realloc(ab, gn, qlen, abpt, simd_p16);
        bits = 16;
    } else {
        simd_p32.inf_min = MAX_OF_THREE(INT32_MIN + abpt->min_mis, INT32_MIN + gap_oe1, INT32_MIN + gap_oe2) + 512 * MAX_OF_TWO(gap_ext1, gap_ext2);
        mem_ret = simd_abpoa_realloc(ab, gn, qlen, abpt, simd_p32);
        bits = 32;
    }
    if (mem_ret) return 0;

    if (bits == 16) {
        if (abpt->gap_mode == ABPOA_LINEAR_GAP) {
            simd_abpoa_lg_align_sequence_to_graph_core(int16_t, simd_p16, SIMDSetOnei16, SIMDMaxi16, \
                    SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16);
        } else if (abpt->gap_mode == ABPOA_AFFINE_GAP) {
            simd_abpoa_ag_align_sequence_to_graph_core(int16_t, simd_p16, SIMDSetOnei16, SIMDMaxi16, \
                    SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, SIMDSetIfEquali16);
        } else if (abpt->gap_mode == ABPOA_CONVEX_GAP) {
            simd_abpoa_cg_align_sequence_to_graph_core(int16_t, simd_p16, SIMDSetOnei16, SIMDMaxi16, \
                    SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, SIMDSetIfEquali16);
        }
    } else { // 2147483647, DP_H/E/F: 32 bits
        if (abpt->gap_mode == ABPOA_LINEAR_GAP) {
            simd_abpoa_lg_align_sequence_to_graph_core(int32_t, simd_p32, SIMDSetOnei32, SIMDMaxi32, \
                    SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32);
        } else if (abpt->gap_mode == ABPOA_AFFINE_GAP) {
            simd_abpoa_ag_align_sequence_to_graph_core(int32_t, simd_p32, SIMDSetOnei32, SIMDMaxi32, \
                    SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, SIMDSetIfEquali32);
        } else if (abpt->gap_mode == ABPOA_CONVEX_GAP) {
            simd_abpoa_cg_align_sequence_to_graph_core(int32_t, simd_p32, SIMDSetOnei32, SIMDMaxi32, \
                    SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, SIMDSetIfEquali32);
        }
    }
#endif
    free(index_map);
    abpt->verbose = old_verbose;
    return 0;
}

#ifndef ABPOA_SIMD_DISPATCH
int simd_abpoa_align_sequence_to_graph(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, abpoa_res_t *res) {
    return simd_abpoa_align_sequence_to_subgraph(ab, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, query, qlen, res);
}
#endif
