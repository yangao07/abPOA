#include <stdio.h>
#include <stdlib.h>
#include "abpoa_graph.h"
#include "abpoa_align.h"
#include "simd_instruction.h"
#include "utils.h"

typedef struct {
    const int reg_n, bits_n, num_of_value, size;
    int inf_min; // based on penalty of mismatch and gap_oe
} SIMD_para_t;

#define SIMDShiftOneNi8  1
#define SIMDShiftOneNi16 2
#define SIMDShiftOneNi32 4
#define SIMDShiftOneNi64 8

#ifdef __AVX512F__
SIMD_para_t _simd_p8  = {512,  8, 64, 64, -1};
SIMD_para_t _simd_p16 = {512, 16, 32, 64, -1};
SIMD_para_t _simd_p32 = {512, 32, 16, 64, -1};
SIMD_para_t _simd_p64 = {512, 64,  8, 64, -1};
#define SIMDTotalBytes 64
#elif defined(__AVX2__)
SIMD_para_t _simd_p8  = {256,  8, 32, 32, -1};
SIMD_para_t _simd_p16 = {256, 16, 16, 32, -1};
SIMD_para_t _simd_p32 = {256, 32,  8, 32, -1};
SIMD_para_t _simd_p64 = {256, 64,  4, 32, -1};
#define SIMDTotalBytes 32
#else
SIMD_para_t _simd_p8  = {128,  8, 16, 16, -1};
SIMD_para_t _simd_p16 = {128, 16,  8, 16, -1};
SIMD_para_t _simd_p32 = {128, 32,  4, 16, -1};
SIMD_para_t _simd_p64 = {128, 64,  2, 16, -1};
#define SIMDTotalBytes 16
#endif

#define print_simd(s, str, score_t) {                               \
    int i; score_t *_a = (score_t*)s;                               \
    printf("%s\t", str);                                            \
    for (i = 0; i < SIMDTotalBytes / (int)sizeof(score_t); ++i) {   \
        printf("%d\t", _a[i]);                                      \
    } printf("\n");                                                 \
}

#define simd_abpoa_print_ag_matrix(score_t) { \
    for (j = 0; j <= target_node_n; ++j) {	\
        printf("index: %d\t", j);	\
        dp_h = DP_HE + j * 2 * dp_sn, dp_e = dp_h + dp_sn;	\
        _dp_h = (int16_t*)dp_h, _dp_e = (int16_t*)dp_e;	\
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {	\
            printf("%d:(%d,%d)\t", i, _dp_h[i], _dp_e[i]);	\
        } printf("\n");	\
    }	\
}

#define simd_abpoa_print_lg_matrix(score_t) { \
    for (j = 0; j <= graph->node_n-2; ++j) {	\
        printf("index: %d\t", j);	\
        dp_h = DP_H + j * dp_sn;	\
        _dp_h = (score_t*)dp_h;	\
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {	\
            printf("%d:(%d)\t", i, _dp_h[i]);	\
        } printf("\n");	\
    }	\
}

/* === workflow of alignment === */
/* a. global:
 * (1) alloc mem 
 * (2) init for first row
 * (3) DP for each row
 * (3.2) if use_ada, update min/max_rank
 * (4) find best_i/j, backtrack
 * b. extend:
 * (1) alloc mem
 * (2) init for first row
 * (3) DP for each row
 * (3.2) find max of current row
 * (3.3) z-drop, set_max_score
 * (3.4) if use_ada, update min/max_rank
 * c. local:
 * (1) alloc mem (zero)
 * (2) init for first row with min_zero
 * (3) DP for each row with min_zero
 * (3.2) find max of current row
 * (3.3) if max > 0 && use_ada, update min/max_rank
 */

#define simd_abpoa_ag_backtrack(score_t, DP_HE, dp_sn, match, mis, gap_e, pre_index, pre_n, backtrack_z, d_sn,                 \
        start_i, start_j, best_i, best_j, qlen, graph, query, n_cigar, graph_cigar) {                                       \
    int i, j, k, pre_i;                                                                                                     \
    SIMDi *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *bz;                                                                          \
    score_t *_dp_h, *_pre_dp_h, *_dp_e, *_pre_dp_e;                                                                         \
    int n_c = 0, s, m_c = 0, id, which, last_which;                                                                         \
    int op_shift[3] = {HOP_OFF_SET, EOP_OFF_SET, FOP_OFF_SET};                                                              \
    score_t d, *_bz; abpoa_cigar_t *cigar = 0;                                                                              \
                                                                                                                            \
    i = best_i, j = best_j, id = abpoa_graph_index_to_node_id(graph, i), last_which = 0;                                    \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, qlen-j, -1, qlen-1);                     \
    while (i > start_i && j > start_j) {                                                                                    \
        bz = backtrack_z + (i-1) * d_sn;                                                                                    \
        _bz = (score_t*)bz;                                                                                                 \
        d = _bz[j];                                                                                                         \
        which = (d >> op_shift[last_which]) & 3;                                                                            \
        if (which == 0) { /* match */                                                                                       \
            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                                            \
            s = graph->node[id].base == query[j-1] ? match : -mis;                                                          \
                                                                                                                            \
            for (k = 0; k < pre_n[i]; ++k) {                                                                                \
                pre_i = pre_index[i][k];                                                                                    \
                dp_h = DP_HE + i * 2 * dp_sn; pre_dp_h = DP_HE + pre_i * 2 * dp_sn;                                         \
                _dp_h = (score_t*)dp_h; _pre_dp_h = (score_t*)pre_dp_h;                                                     \
                if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                                       \
                    i = pre_i;                                                                                              \
                    break;                                                                                                  \
                }                                                                                                           \
            }                                                                                                               \
            id = abpoa_graph_index_to_node_id(graph, i);                                                                    \
            j--; last_which = which;                                                                                        \
        } else if (which == 1) { /* deletion */                                                                             \
            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                                              \
            if (last_which == 1) {                                                                                          \
                for (k = 0; k < pre_n[i]; ++k) {                                                                            \
                    pre_i = pre_index[i][k];                                                                                \
                    dp_e = DP_HE + (i * 2 + 1) * dp_sn; pre_dp_e = DP_HE + (pre_i * 2 + 1) * dp_sn;                         \
                    _dp_e = (score_t*)dp_e; _pre_dp_e = (score_t*)pre_dp_e;                                                 \
                    if (_pre_dp_e[j] - gap_e == _dp_e[j]) {                                                                 \
                        i = pre_i;                                                                                          \
                        break;                                                                                              \
                    }                                                                                                       \
                }                                                                                                           \
            } else if (last_which == 0) {                                                                                   \
                for (k = 0; k < pre_n[i]; ++k) {                                                                            \
                    pre_i = pre_index[i][k];                                                                                \
                    dp_h = DP_HE + i * 2 * dp_sn; pre_dp_e = DP_HE + (pre_i * 2 + 1) * dp_sn;                               \
                    _dp_h = (score_t*)dp_h; _pre_dp_e = (score_t*)pre_dp_e;                                                 \
                    if (_pre_dp_e[j] == _dp_h[j]) {                                                                         \
                        i = pre_i;                                                                                          \
                        break;                                                                                              \
                    }                                                                                                       \
                }                                                                                                           \
            } else {                                                                                                        \
                _err_fatal_simple(__func__, "\nunexpected cigar op.\n");                                                    \
            }                                                                                                               \
            id = abpoa_graph_index_to_node_id(graph, i);                                                                    \
            last_which = which;                                                                                             \
        } else { /* insertion */                                                                                            \
            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1);                                              \
            j--; last_which = which;                                                                                        \
        }                                                                                                                   \
    }                                                                                                                       \
    if (j > start_j) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, j-start_j, -1, j-1);                       \
    /* reverse cigar */                                                                                                     \
    *graph_cigar = abpoa_reverse_cigar(n_c, cigar);                                                                         \
    *n_cigar = n_c;                                                                                                         \
    /* DEBUG */                                                                                                             \
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/                                                                        \
}

#define simd_abpoa_lg_backtrack(score_t, DP_H, dp_sn, match, mis, gap_e, pre_index, pre_n,                                  \
        start_i, start_j, best_i, best_j, qlen, graph, query, n_cigar, graph_cigar) {                                       \
    int i, j, k, pre_i;                                                                                                     \
    SIMDi *dp_h, *pre_dp_h; score_t *_dp_h, *_pre_dp_h;                                                                     \
    int n_c = 0, s, m_c = 0, hit, id;                                                                                       \
    abpoa_cigar_t *cigar = 0;                                                                                               \
                                                                                                                            \
    i = best_i, j = best_j, id = abpoa_graph_index_to_node_id(graph, i);                                                    \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, qlen-j, -1, qlen-1);                     \
    while (i > start_i && j > start_j) {                                                                                    \
        hit = 0;                                                                                                            \
        s = graph->node[id].base == query[j-1] ? match : -mis;                                                              \
        for (k = 0; k < pre_n[i]; ++k) {                                                                                    \
            pre_i = pre_index[i][k];                                                                                        \
            dp_h = DP_H + i * dp_sn; pre_dp_h = DP_H + pre_i * dp_sn;                                                       \
            _dp_h = (score_t*)dp_h; _pre_dp_h = (score_t*)pre_dp_h;                                                         \
            if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                                           \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                                        \
                i = pre_i;                                                                                                  \
                id = abpoa_graph_index_to_node_id(graph, i);                                                                \
                j--;                                                                                                        \
                hit = 1;                                                                                                    \
                break;                                                                                                      \
            } else if (_pre_dp_h[j] - gap_e == _dp_h[j]) {                                                                  \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                                          \
                i = pre_i;                                                                                                  \
                id = abpoa_graph_index_to_node_id(graph, i);                                                                \
                hit = 1;                                                                                                    \
                break;                                                                                                      \
            }                                                                                                               \
        }                                                                                                                   \
        if (hit == 0 && _dp_h[j-1] - gap_e == _dp_h[j]) {                                                                   \
            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1);                                              \
            j--;                                                                                                            \
            continue;                                                                                                       \
        }                                                                                                                   \
    }                                                                                                                       \
    if (j > start_j) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, j-start_j, -1, j-1);                       \
    /* reverse cigar */                                                                                                     \
    *graph_cigar = abpoa_reverse_cigar(n_c, cigar);                                                                         \
    *n_cigar = n_c;                                                                                                         \
    /* DEBUG */                                                                                                             \
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/                                                                        \
}

#define simd_abpoa_ag_local_backtrack(score_t, DP_HE, dp_sn, match, mis, gap_e, pre_index, pre_n, backtrack_z, d_sn,        \
        start_i, start_j, best_i, best_j, qlen, graph, query, n_cigar, graph_cigar) {                                       \
    int i, j, k, pre_i;                                                                                                     \
    SIMDi *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *bz;                                                                          \
    score_t *_dp_h, *_pre_dp_h, *_dp_e, *_pre_dp_e;                                                                         \
    int n_c = 0, s, m_c = 0, id, which, last_which;                                                                         \
    int op_shift[3] = {HOP_OFF_SET, EOP_OFF_SET, FOP_OFF_SET};                                                              \
    score_t d, *_bz; abpoa_cigar_t *cigar = 0;                                                                              \
                                                                                                                            \
    i = best_i, j = best_j, id = abpoa_graph_index_to_node_id(graph, i), last_which = 0;                                    \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, qlen-j, -1, qlen-1);                     \
    while (1) {                                                                                                             \
        dp_h = DP_HE + i * 2 * dp_sn; _dp_h = (score_t*)dp_h;                                                               \
        if (_dp_h[j] == 0) break;                                                                                           \
        bz = backtrack_z + (i-1) * d_sn;                                                                                    \
        _bz = (score_t*)bz;                                                                                                 \
        d = _bz[j];                                                                                                         \
        which = (d >> op_shift[last_which]) & 3;                                                                            \
        if (which == 0) { /* match */                                                                                       \
            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                                            \
            s = graph->node[id].base == query[j-1] ? match : -mis;                                                          \
                                                                                                                            \
            for (k = 0; k < pre_n[i]; ++k) {                                                                                \
                pre_i = pre_index[i][k];                                                                                    \
                pre_dp_h = DP_HE + pre_i * 2 * dp_sn;                                                                       \
                _pre_dp_h = (score_t*)pre_dp_h;                                                                             \
                if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                                       \
                    i = pre_i;                                                                                              \
                    break;                                                                                                  \
                }                                                                                                           \
            }                                                                                                               \
            id = abpoa_graph_index_to_node_id(graph, i);                                                                    \
            j--; last_which = which;                                                                                        \
        } else if (which == 1) { /* deletion */                                                                             \
            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                                              \
            if (last_which == 1) {                                                                                          \
                for (k = 0; k < pre_n[i]; ++k) {                                                                            \
                    pre_i = pre_index[i][k];                                                                                \
                    dp_e = DP_HE + (i * 2 + 1) * dp_sn; pre_dp_e = DP_HE + (pre_i * 2 + 1) * dp_sn;                         \
                    _dp_e = (score_t*)dp_e; _pre_dp_e = (score_t*)pre_dp_e;                                                 \
                    if (_pre_dp_e[j] - gap_e == _dp_e[j]) {                                                                 \
                        i = pre_i;                                                                                          \
                        break;                                                                                              \
                    }                                                                                                       \
                }                                                                                                           \
            } else if (last_which == 0) {                                                                                   \
                for (k = 0; k < pre_n[i]; ++k) {                                                                            \
                    pre_i = pre_index[i][k];                                                                                \
                    pre_dp_e = DP_HE + (pre_i * 2 + 1) * dp_sn;                                                             \
                    _pre_dp_e = (score_t*)pre_dp_e;                                                                         \
                    if (_pre_dp_e[j] == _dp_h[j]) {                                                                         \
                        i = pre_i;                                                                                          \
                        break;                                                                                              \
                    }                                                                                                       \
                }                                                                                                           \
            } else {                                                                                                        \
                _err_fatal_simple(__func__, "\nunexpected cigar op.\n");                                                    \
            }                                                                                                               \
            id = abpoa_graph_index_to_node_id(graph, i);                                                                    \
            last_which = which;                                                                                             \
        } else { /* insertion */                                                                                            \
            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1);                                              \
            j--; last_which = which;                                                                                        \
        }                                                                                                                   \
    }                                                                                                                       \
    if (j > start_j) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, j-start_j, -1, j-1);                       \
    /* reverse cigar */                                                                                                     \
    *graph_cigar = abpoa_reverse_cigar(n_c, cigar);                                                                         \
    *n_cigar = n_c;                                                                                                         \
    /* DEBUG */                                                                                                             \
    abpoa_print_cigar(n_c, *graph_cigar, graph);                                                                            \
}

#define simd_abpoa_lg_local_backtrack(score_t, DP_H, dp_sn, match, mis, gap_e, pre_index, pre_n,                            \
        start_i, start_j, best_i, best_j, qlen, graph, query, n_cigar, graph_cigar) {                                       \
    int i, j, k, pre_i;                                                                                                     \
    SIMDi *dp_h, *pre_dp_h; score_t *_dp_h, *_pre_dp_h;                                                                     \
    int n_c = 0, s, m_c = 0, hit, id;                                                                                       \
    abpoa_cigar_t *cigar = 0;                                                                                               \
                                                                                                                            \
    i = best_i, j = best_j, id = abpoa_graph_index_to_node_id(graph, i);                                                    \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, qlen-j, -1, qlen-1);                     \
    while (1) {                                                                                                             \
        dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h;                                                                    \
        if (_dp_h[j] == 0) break;                                                                                           \
        s = graph->node[id].base == query[j-1] ? match : -mis;                                                              \
        hit = 0;                                                                                                            \
        for (k = 0; k < pre_n[i]; ++k) {                                                                                    \
            pre_i = pre_index[i][k];                                                                                        \
            dp_h = DP_H + i * dp_sn; pre_dp_h = DP_H + pre_i * dp_sn;                                                       \
            _dp_h = (score_t*)dp_h; _pre_dp_h = (score_t*)pre_dp_h;                                                         \
            if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                                           \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                                        \
                i = pre_i;                                                                                                  \
                id = abpoa_graph_index_to_node_id(graph, i);                                                                \
                j--;                                                                                                        \
                hit = 1;                                                                                                    \
                break;                                                                                                      \
            } else if (_pre_dp_h[j] - gap_e == _dp_h[j]) {                                                                  \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                                          \
                i = pre_i;                                                                                                  \
                id = abpoa_graph_index_to_node_id(graph, i);                                                                \
                hit = 1;                                                                                                    \
                break;                                                                                                      \
            }                                                                                                               \
        }                                                                                                                   \
        if (hit == 0 && _dp_h[j-1] - gap_e == _dp_h[j]) {                                                                   \
            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1);                                              \
            j--;                                                                                                            \
            continue;                                                                                                       \
        }                                                                                                                   \
    }                                                                                                                       \
    if (j > start_j) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CSOFT_CLIP, j-start_j, -1, j-1);                       \
    /* reverse cigar */                                                                                                     \
    *graph_cigar = abpoa_reverse_cigar(n_c, cigar);                                                                         \
    *n_cigar = n_c;                                                                                                         \
    /* DEBUG */                                                                                                             \
    abpoa_print_cigar(n_c, *graph_cigar, graph);                                                                            \
}

#define simd_abpoa_ag_var(score_t, sp, SIMDSetOne) \
    abpoa_graph_t *graph = ab->abg; abpoa_simd_matrix_t *abm = ab->abm;           \
    int **pre_index, *pre_n, pre_i; \
    int matrix_row_n = graph->node_n, matrix_col_n = qlen + 1;   \
    int i, j, k, *dp_beg, *dp_beg_sn, *dp_end, *dp_end_sn, node_id, index_i, q_i;    \
    int pn , qp_sn,d_sn, dp_sn, size; /* pn: value per SIMDi, qp_sn/dp_sn/d_sn: segmented length*/  \
    int beg, end, beg_sn, end_sn, _beg_sn, _end_sn, pre_beg_sn, pre_beg, pre_end, pre_end_sn, sn_i; \
    \
    SIMDi *qp, *qi;                                                        \
    SIMDi *DP_HE, *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *dp_f, tmp_e;                                 \
    SIMDi *backtrack_z, *z; /* backtrack cell: f<<4|e<<2|h, MATCH:0, DELETION:1, INSERTION:2 */     \
    SIMDi *hd, fd, ed, hm, he, hf, ee, em, ff, fm; uint8_t m0=0x0, e1=0x1, f2=0x2;                  \
    \
    score_t *_dp_h, *_qi, *_dp_e, *_dp_f, best_score = sp.inf_min, inf_min = sp.inf_min;                  \
    int *mat = abpt->mat, best_i = 0, best_j = 0, max, max_i;                                       \
    score_t gap_open = abpt->gap_open, gap_ext = abpt->gap_ext, gap_OE = abpt->gap_open + abpt->gap_ext;    \
    \
    int w = abpt->bw < 0 ? qlen : abpt->bw; /* when w < 0, do whole global */   \
    SIMDi zero = SIMDSetZeroi(), min_inf = SIMDSetOne(inf_min), gap_e = SIMDSetOne(gap_ext), gap_oe = SIMDSetOne(gap_OE); \
\
{   /* allocate memory */   \
    pn = sp.num_of_value, size = sp.size; qp_sn = d_sn = dp_sn = (matrix_col_n + pn - 1) / pn;    \
    qp = abm->s_mem;   \
    /* generate the query profile */    \
    for (i = 0; i < qp_sn * abpt->m; ++i) qp[i] = min_inf;  \
    for (k = 0; k < abpt->m; ++k) { /* SIMD parallelization */ \
        int *p = &mat[k * abpt->m]; \
        score_t *_qp = (score_t*)(qp + k * qp_sn);  \
        _qp[0] = 0;     \
        for (j = 0; j < qlen; ++j) _qp[j+1] = (score_t)p[query[j]]; \
        for (j = qlen+1; j < qp_sn * pn; ++j) _qp[j] = 0;   \
    }   \
    DP_HE = qp + qp_sn * abpt->m; /* H|E|H|E */ \
    dp_f = DP_HE + 2 * dp_sn * matrix_row_n;    \
    if (abpt->use_ada || abpt->align_mode == ABPOA_LOCAL_MODE || abpt->align_mode == ABPOA_EXTEND_MODE) { /* query index */   \
        qi = dp_f + dp_sn;  \
        _qi = (score_t*)qi;    \
        for (i = 0; i <= qlen; ++i) _qi[i] = i; \
        for (i = qlen+1; i < (qlen/pn+1) * pn; ++i) _qi[i] = -1;    \
    }   \
    if (abpt->ret_cigar) hd = qi + dp_sn, backtrack_z = hd + d_sn;  \
    else if (abpt->use_ada) hd = qi + dp_sn;    \
    /* for backtrack */ \
    if (abpt->ret_cigar || abpt->use_ada) {   \
        hm = SIMDSetOne(m0 << HOP_OFF_SET), he = SIMDSetOne(e1 << HOP_OFF_SET), hf = SIMDSetOne(f2 << HOP_OFF_SET);    \
        em = SIMDSetOne(m0 << EOP_OFF_SET), ee = SIMDSetOne(e1 << EOP_OFF_SET), fm = SIMDSetOne(m0 << FOP_OFF_SET), ff = SIMDSetOne(f2 << FOP_OFF_SET); \
    }   \
    dp_beg = abm->dp_beg, dp_end = abm->dp_end, dp_beg_sn = abm->dp_beg_sn, dp_end_sn = abm->dp_end_sn; \
    if (abpt->bw >= 0) { \
        for (i = 1; i < graph->node_n; ++i) {   \
            node_id = graph->index_to_node_id[i];   \
            graph->node_id_to_min_rank[node_id] = graph->node_n;    \
            graph->node_id_to_max_rank[node_id] = 0;    \
        }   \
    }   \
    /* index of pre-node */ \
    pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));   \
    pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));    \
    for (i = 0; i < graph->node_n; ++i) {   \
        node_id = abpoa_graph_index_to_node_id(graph, i); /* i: node index */ \
        pre_n[i] = graph->node[node_id].in_edge_n;  \
        pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));   \
        for (j = 0; j < pre_n[i]; ++j) {    \
            pre_index[i][j] = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);   \
        }   \
    }   \
}

#define simd_abpoa_free_var { \
    for (i = 0; i < graph->node_n; ++i) free(pre_index[i]); free(pre_index); free(pre_n);	\
}	\

#define simd_abpoa_ag_first_row(score_t) { \
    /* fill the first row */    \
    if (abpt->bw >= 0) {    \
        graph->node_id_to_min_rank[0] = graph->node_id_to_max_rank[0] = 0;  \
        for (i = 0; i < graph->node[0].out_edge_n; ++i) { /* set min/max rank for next_id */ \
            int out_id = graph->node[0].out_id[i];  \
            graph->node_id_to_min_rank[out_id] = graph->node_id_to_max_rank[out_id] = 1;    \
        }   \
        dp_beg[0] = GET_DP_BEGIN(graph, w, 0, qlen), dp_end[0] = GET_DP_END(graph, w, 0, qlen);    \
    } else { \
        dp_beg[0] = 0, dp_end[0] = qlen;    \
    }   \
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;   \
    dp_h = DP_HE, dp_e = dp_h + dp_sn;  \
    _dp_h = (score_t*)dp_h, _dp_e = (score_t*)dp_e; \
    _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);  \
}

#define simd_abpoa_ag_first_dp(score_t) { \
    simd_abpoa_ag_first_row(score_t);   \
    for (i = 0; i <= _end_sn; ++i) {    \
        dp_h[i] = min_inf; dp_e[i] = min_inf;   \
    }   \
    for (i = 1; i <= dp_end[0]; ++i) { /* no SIMD parallelization */\
        _dp_h[i] = -gap_open - gap_ext * i; \
    } _dp_h[0] = 0; _dp_e[0] = -(gap_OE);   \
}

#define simd_abpoa_ag_local_first_dp(score_t) { \
    simd_abpoa_ag_first_row(score_t);    \
    for (i = 0; i <= _end_sn; ++i) {    \
        dp_h[i] = zero; dp_e[i] = zero; \
    }   \
}

#define simd_abpoa_ag_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater) { \
    node_id = abpoa_graph_index_to_node_id(graph, index_i); \
    SIMDi *q = qp + graph->node[node_id].base * qp_sn;  \
    dp_h = DP_HE + index_i * 2 * dp_sn, dp_e = dp_h + dp_sn;    \
    _dp_h = (score_t*)dp_h, _dp_e = (score_t*)dp_e, _dp_f = (score_t*)dp_f; \
    if (abpt->bw < 0) beg = dp_beg[index_i] = 0, end = dp_end[index_i] = qlen;  \
    else beg = dp_beg[index_i] = GET_DP_BEGIN(graph, w, node_id, qlen), end = dp_end[index_i] = GET_DP_END(graph, w, node_id, qlen);    \
    beg_sn = dp_beg_sn[index_i] = (dp_beg[index_i])/pn; end_sn = dp_end_sn[index_i] = (dp_end[index_i])/pn; \
    /* loop query */    \
    /* init h, e */     \
    _beg_sn = MAX_OF_TWO(beg_sn-1, 0);  \
    _end_sn = MIN_OF_TWO(end_sn+1, dp_sn-1);    \
    for (i = _beg_sn; i <= _end_sn; ++i) { /* SIMD parallelization */   \
        dp_h[i] = min_inf; dp_e[i] = min_inf; dp_f[i] = min_inf;    \
    }   \
    /* get max m and e */   \
    for (i = 0; i < pre_n[index_i]; ++i) {  \
        pre_i = pre_index[index_i][i];  \
        pre_dp_h = DP_HE + pre_i * 2 * dp_sn, pre_dp_e = pre_dp_h + dp_sn;  \
        pre_beg = dp_beg[pre_i]; pre_end = dp_end[pre_i];   \
        pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];   \
        /* set M from (pre_i, q_i-1) */ \
        _beg_sn = MAX_OF_THREE(0, pre_beg_sn, (beg-1)/pn);  \
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);    \
        SIMDi first = SIMDShiftRight(min_inf, SIMDTotalBytes-SIMDShiftOneN); \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */ \
            SIMDi remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);  \
            dp_h[sn_i] = SIMDMax(SIMDOri(first, remain), dp_h[sn_i]); \
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN); \
        }   \
        /* set E from (pre_i, q_i) */   \
        _beg_sn = MAX_OF_TWO(pre_beg_sn, beg_sn); _end_sn = MIN_OF_TWO(pre_end_sn, end_sn); \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */  \
            dp_e[sn_i] = SIMDMax(pre_dp_e[sn_i], dp_e[sn_i]);  \
        }   \
    }   \
    if (abpt->ret_cigar) {   \
        z = backtrack_z + (index_i-1) * d_sn;   \
        /* compare M, E, and F */   \
        for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */\
            SIMDGetIfGreater(hd[sn_i], dp_h[sn_i], dp_e[sn_i], SIMDAdd(dp_h[sn_i], q[sn_i]), he, hm);    \
        }   \
        /* set F from (index_i, q_i-1) */   \
        for (q_i = beg+1; q_i <= end; ++q_i) { /* no SIMD parallelization */    \
            _dp_f[q_i] = MAX_OF_TWO(_dp_h[q_i-1] - gap_OE, _dp_f[q_i-1] - gap_ext); \
        }   \
        for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */\
            /* h, hd for current cell */\
            SIMDGetIfGreater(hd[sn_i], dp_h[sn_i], dp_f[sn_i], dp_h[sn_i], hf, hd[sn_i]); \
            /* e, ed for next cell */   \
            tmp_e = SIMDSub(dp_h[sn_i], gap_oe);  \
            SIMDGetIfGreater(ed, dp_e[sn_i], SIMDSub(dp_e[sn_i], gap_e), tmp_e, ee, em);   \
            /* fd for next cell */  \
            SIMDSetIfGreater(fd, SIMDSub(dp_f[sn_i], gap_e), tmp_e, ff, fm);  \
            z[sn_i] = SIMDOri(SIMDOri(hd[sn_i], ed), fd);  \
        }   \
    } else {    \
        /* compare M, E, and F */   \
        for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */\
            if (abpt->use_ada) {    \
                SIMDGetIfGreater(hd[sn_i], dp_h[sn_i], dp_e[sn_i], SIMDAdd(dp_h[sn_i], q[sn_i]), he, hm);    \
            } else {    \
                dp_h[sn_i] = SIMDMax(SIMDAdd(dp_h[sn_i], q[sn_i]), dp_e[sn_i]);  \
            }   \
        }   \
        /* set F from (index_i, q_i-1) */   \
        _dp_f[beg] = inf_min;   \
        for (q_i = beg+1; q_i <= end; ++q_i) { /* no SIMD parallelization */\
            _dp_f[q_i] = MAX_OF_TWO(_dp_h[q_i-1] - gap_OE, _dp_f[q_i-1] - gap_ext); \
        }   \
        for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */    \
            /* h for current cell */ \
            dp_h[sn_i] = SIMDMax(dp_h[sn_i], dp_f[sn_i]); \
            /* e for next cell */   \
            dp_e[sn_i] = SIMDMax(SIMDSub(dp_e[sn_i], gap_e), SIMDSub(dp_h[sn_i], gap_oe));  \
        }   \
    } \
}

#define simd_abpoa_ag_local_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater) { \
    node_id = abpoa_graph_index_to_node_id(graph, index_i); \
    SIMDi *q = qp + graph->node[node_id].base * qp_sn;  \
    dp_h = DP_HE + index_i * 2 * dp_sn, dp_e = dp_h + dp_sn;    \
    _dp_h = (score_t*)dp_h, _dp_e = (score_t*)dp_e, _dp_f = (score_t*)dp_f; \
    if (abpt->bw < 0) beg = dp_beg[index_i] = 0, end = dp_end[index_i] = qlen;  \
    else beg = dp_beg[index_i] = GET_DP_BEGIN(graph, w, node_id, qlen), end = dp_end[index_i] = GET_DP_END(graph, w, node_id, qlen);    \
    beg_sn = dp_beg_sn[index_i] = (dp_beg[index_i])/pn; end_sn = dp_end_sn[index_i] = (dp_end[index_i])/pn; \
    /* loop query */    \
    /* init h, e */ \
    _beg_sn = MAX_OF_TWO(beg_sn-1, 0);  \
    _end_sn = MIN_OF_TWO(end_sn+1, dp_sn-1);\
    for (i = _beg_sn; i <= _end_sn; ++i) { /* SIMD parallelization */ \
        dp_h[i] = zero; dp_e[i] = zero; dp_f[i] = zero;\
    }\
    /* get max m and e */ \
    for (i = 0; i < pre_n[index_i]; ++i) {\
        pre_i = pre_index[index_i][i];\
        pre_dp_h = DP_HE + pre_i * 2 * dp_sn, pre_dp_e = pre_dp_h + dp_sn;\
        pre_beg = dp_beg[pre_i]; pre_end = dp_end[pre_i];\
        pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];\
        /* set M from (pre_i, q_i-1) */ \
        _beg_sn = MAX_OF_THREE(0, pre_beg_sn, (beg-1)/pn);\
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);\
        SIMDi first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);\
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */\
            SIMDi remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);\
            dp_h[sn_i] = SIMDMax(SIMDOri(first, remain), dp_h[sn_i]); \
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN); \
        }   \
        /* set E from (pre_i, q_i) */ \
        _beg_sn = MAX_OF_TWO(pre_beg_sn, beg_sn); _end_sn = MIN_OF_TWO(pre_end_sn, end_sn); \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */ \
            dp_e[sn_i] = SIMDMax(pre_dp_e[sn_i], dp_e[sn_i]);  \
        }   \
    }   \
    if (abpt->ret_cigar) {\
        z = backtrack_z + (index_i-1) * d_sn;\
        /* compare M, E, and F */ \
        for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */ \
            SIMDGetIfGreater(hd[sn_i], dp_h[sn_i], dp_e[sn_i], SIMDAdd(dp_h[sn_i], q[sn_i]), he, hm);\
        }\
        /* set F from (index_i, q_i-1) */ \
        for (q_i = beg+1; q_i <= end; ++q_i) { /* no SIMD parallelization */ \
            _dp_f[q_i] = MAX_OF_THREE(0, _dp_h[q_i-1] - gap_OE, _dp_f[q_i-1] - gap_ext);\
        }                                       \
        for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */ \
            /* h, hd for current cell */ \
            SIMDGetIfGreater(hd[sn_i], dp_h[sn_i], dp_f[sn_i], dp_h[sn_i], hf, hd[sn_i]);\
            /* e, ed for next cell */ \
            tmp_e = SIMDSub(dp_h[sn_i], gap_oe);\
            SIMDGetIfGreater(ed, dp_e[sn_i], SIMDSub(dp_e[sn_i], gap_e), tmp_e, ee, em);\
            dp_e[sn_i] = SIMDMax(dp_e[sn_i], zero); \
            /* fd for next cell */\
            SIMDSetIfGreater(fd, SIMDSub(dp_f[sn_i], gap_e), tmp_e, ff, fm);\
            z[sn_i] = SIMDOri(SIMDOri(hd[sn_i], ed), fd);\
        }\
    } else {\
        /* compare M, E, and F */ \
        for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */ \
            if (abpt->use_ada) {\
                SIMDGetIfGreater(hd[sn_i], dp_h[sn_i], dp_e[sn_i], SIMDAdd(dp_h[sn_i], q[sn_i]), he, hm);\
            } else {\
                dp_h[sn_i] = SIMDMax(SIMDAdd(dp_h[sn_i], q[sn_i]), dp_e[sn_i]);\
            }\
        }\
        /* set F from (index_i, q_i-1) */ \
        _dp_f[beg] = 0;\
        for (q_i = beg+1; q_i <= end; ++q_i) { /* no SIMD parallelization */ \
            _dp_f[q_i] = MAX_OF_THREE(0, _dp_h[q_i-1] - gap_OE, _dp_f[q_i-1] - gap_ext);\
        }\
        for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */ \
            /* h for current cell */\
            dp_h[sn_i] = SIMDMax(dp_h[sn_i], dp_f[sn_i]);\
            /* e for next cell */ \
            dp_e[sn_i] = SIMDMax(zero, SIMDMax(SIMDSub(dp_e[sn_i], gap_e), SIMDSub(dp_h[sn_i], gap_oe)));\
        }\
    }\
}

#define simd_abpoa_max(score_t, SIMDSetIfGreater, SIMDGetIfGreater) { \
    /* select max dp_h */   \
    max = inf_min, max_i = -1;   \
    SIMDi a = dp_h[end_sn], b = qi[end_sn];   \
    if (end_sn == qlen / pn) SIMDSetIfGreater(a, zero, b, min_inf, a);  \
    for (i = beg_sn; i < end_sn; ++i) { \
        SIMDGetIfGreater(b, a, dp_h[i], a, qi[i], b);   \
    }   \
    _dp_h = (score_t*)&a, _qi = (score_t*)&b; \
    for (i = 0; i < pn; ++i) {  \
        if (_dp_h[i] > max) {   \
            max = _dp_h[i]; \
            max_i = _qi[i];  \
        }   \
    }   \
}

#define simd_abpoa_ag_global_get_max(score_t) {	\
    int in_id, in_index;	\
    for (i = 0; i < graph->node[ABPOA_SINK_NODE_ID].in_edge_n; ++i) { 	\
        in_id = graph->node[ABPOA_SINK_NODE_ID].in_id[i];	\
        in_index = abpoa_graph_node_id_to_index(graph, in_id);	\
        dp_h = DP_HE + in_index * 2 * dp_sn;	\
        _dp_h = (score_t*)dp_h;	\
        _set_max_score(best_score, best_i, best_j, _dp_h[qlen], in_index, qlen);	\
    }	\
}

#define simd_abpoa_ada_rank_from_maxpre { \
    /* set min/max_rank for next nodes */ \
    graph->node_id_to_max_rank[node_id] = graph->node_id_to_max_rank[max_pre_id] + 1;\
    graph->node_id_to_min_rank[node_id] = graph->node_id_to_min_rank[max_pre_id] + 1;\
    int max_out_rank = graph->node_id_to_max_rank[max_pre_id] + 2;\
    int min_out_rank = graph->node_id_to_min_rank[max_pre_id] + 2;\
    for (i = 0; i < graph->node[node_id].out_edge_n; ++i) {\
        int out_node_id = graph->node[node_id].out_id[i];\
        if (max_out_rank > graph->node_id_to_max_rank[out_node_id]) graph->node_id_to_max_rank[out_node_id] = max_out_rank;\
        if (min_out_rank < graph->node_id_to_min_rank[out_node_id]) graph->node_id_to_min_rank[out_node_id] = min_out_rank;\
    }   \
} 

#define simd_abpoa_rank { \
    /* set min/max_rank for next nodes */\
    int max_out_rank = graph->node_id_to_max_rank[node_id] + 1;\
    int min_out_rank = graph->node_id_to_min_rank[node_id] + 1;\
    for (i = 0; i < graph->node[node_id].out_edge_n; ++i) {\
        int out_node_id = graph->node[node_id].out_id[i];\
        if (max_out_rank > graph->node_id_to_max_rank[out_node_id]) graph->node_id_to_max_rank[out_node_id] = max_out_rank;\
        if (min_out_rank < graph->node_id_to_min_rank[out_node_id]) graph->node_id_to_min_rank[out_node_id] = min_out_rank;\
    }\
}

#define simd_abpoa_ag_ada_rank(score_t) { \
    /* determine max pre_i, then determin current rank */   \
    int which, s, max_pre_i=-1, max_pre_id; \
    score_t *_hd = (score_t*)hd, *_pre_dp_h, *_pre_dp_e;    \
    which = (_hd[max_i] >> HOP_OFF_SET) & 3;    \
    if (which == 0) { /* match */   \
        s = graph->node[node_id].base == query[max_i-1] ? abpt->match : -abpt->mismatch;    \
        for (k = 0; k < pre_n[index_i]; ++k) {  \
            pre_i = pre_index[index_i][k];  \
            dp_h = DP_HE + index_i * 2 * dp_sn; pre_dp_h = DP_HE + pre_i * 2 * dp_sn;   \
            _dp_h = (score_t*)dp_h; _pre_dp_h = (score_t*)pre_dp_h; \
            if (_pre_dp_h[max_i-1] + s == _dp_h[max_i]) {   \
                max_pre_i = pre_i;  \
                break;  \
            }   \
        }   \
        max_pre_id = abpoa_graph_index_to_node_id(graph, max_pre_i);\
    } else if (which == 1) { /* deletion */\
        for (k = 0; k < pre_n[index_i]; ++k) {\
            pre_i = pre_index[index_i][k];\
            dp_h = DP_HE + index_i * 2 * dp_sn; pre_dp_e = DP_HE + (pre_i * 2 + 1) * dp_sn;\
            _dp_h = (score_t*)dp_h; _pre_dp_e = (score_t*)pre_dp_e;\
            if (_pre_dp_e[max_i] == _dp_h[max_i]) {\
                max_pre_i = pre_i;\
                break;\
            } \
        }\
        max_pre_id = abpoa_graph_index_to_node_id(graph, max_pre_i);\
    } else { /* insertion */ \
        err_fatal_simple("Unexpected cigar op.\n");\
    }  \
    simd_abpoa_ada_rank_from_maxpre;    \
}

#define simd_abpoa_ag_get_cigar(score_t) {  \
    if (abpt->ret_cigar) {   \
        simd_abpoa_ag_backtrack(score_t, DP_HE, dp_sn, (score_t)abpt->match, (score_t)abpt->mismatch, (score_t)abpt->gap_ext, pre_index, pre_n, backtrack_z, d_sn, 0, 0, best_i, best_j, qlen, graph, query, n_cigar, graph_cigar); \
    }   \
}

#define simd_abpoa_ag_local_get_cigar(score_t) {    \
    if (abpt->ret_cigar) {   \
        simd_abpoa_ag_local_backtrack(score_t, DP_HE, dp_sn, (score_t)abpt->match, (score_t)abpt->mismatch, (score_t)abpt->gap_ext, pre_index, pre_n, backtrack_z, d_sn, 0, 0, best_i, best_j, qlen, graph, query, n_cigar, graph_cigar);   \
    }   \
}

// affine gap penalty: gap_open > 0
#define simd_abpoa_ag_global_align_sequence_with_graph_core(score_t, ab, query, qlen, abpt, n_cigar, graph_cigar, sp, \
    SIMDSetOne, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, b_score) { \
    simd_abpoa_ag_var(score_t, sp, SIMDSetOne); \
    simd_abpoa_ag_first_dp(score_t);    \
    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {    \
        simd_abpoa_ag_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater);  \
        if (abpt->bw >= 0) {    \
            if (abpt->use_ada) { /* set current rank, set min/max rank for next nodes */   \
                simd_abpoa_max(score_t, SIMDSetIfGreater, SIMDGetIfGreater);    \
                /*DEBUG*/   \
                /*printf("adaMax: %d, %d\n", max_i, max);*/ \
                simd_abpoa_ag_ada_rank(score_t);    \
            } else {    \
                simd_abpoa_rank;    \
            }   \
        }   \
    }   \
    simd_abpoa_ag_global_get_max(score_t);  \
    /*DEBUG*/   \
    /* simd_abpoa_print_ag_matrix(score_t);
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score); */ \
    simd_abpoa_ag_get_cigar(score_t);   \
    simd_abpoa_free_var; \
    b_score = best_score;  \
} 

#define simd_abpoa_ag_extend_align_sequence_with_graph_core(score_t, graph, query, qlen, abpt, n_cigar, graph_cigar, sp,   \
    SIMDSetOne, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, b_score) { \
    simd_abpoa_ag_var(score_t, sp, SIMDSetOne); \
    simd_abpoa_ag_first_dp(score_t);    \
    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {    \
        simd_abpoa_ag_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater);  \
        simd_abpoa_max(score_t, SIMDSetIfGreater, SIMDGetIfGreater);    \
        /* DEBUG*/  \
        /* printf("extendMax: %d, %d\n", max_i, max); */    \
        /* apply z-drop strategy */ \
        if (abpt->zdrop > 0 && best_score - max > abpt->zdrop + gap_ext * abs((best_i-index_i)-(best_j-max_i))) \
            break;  \
        _set_max_score(best_score, best_i, best_j, max, index_i, max_i);    \
        /* set adaptive band */     \
        if (abpt->bw >= 0) { \
            if (abpt->use_ada) { /* set current rank, set min/max rank for next nodes */  \
                simd_abpoa_ag_ada_rank(score_t);    \
            } else {    \
                simd_abpoa_rank;    \
            }   \
        }   \
    }   \
    /*DEBUG*/   \
    /*simd_abpoa_print_ag_matrix(score_t);
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);*/   \
    simd_abpoa_ag_get_cigar(score_t);   \
    simd_abpoa_free_var; \
    b_score = best_score;  \
}

#define simd_abpoa_ag_local_align_sequence_with_graph_core(score_t, graph, query, qlen, abpt, n_cigar, graph_cigar, sp,    \
    SIMDSetOne, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, b_score) { \
    simd_abpoa_ag_var(score_t, sp, SIMDSetOne);     \
    gap_open = gap_open+0;/*XXX*/   \
    simd_abpoa_ag_local_first_dp(score_t);  \
    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {    \
        simd_abpoa_ag_local_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater);    \
        simd_abpoa_max(score_t, SIMDSetIfGreater, SIMDGetIfGreater);    \
        /*DEBUG*/   \
        /*printf("localMax: %d, %d\n", max_i, max);*/   \
        _set_max_score(best_score, best_i, best_j, max, index_i, max_i);    \
        /* set adaptive band */     \
        if (abpt->bw >= 0) { \
            if (max > 0 && abpt->use_ada) { /* set current rank, set min/max rank for next nodes */    \
                simd_abpoa_ag_ada_rank(score_t);    \
            } else {    \
                simd_abpoa_rank;    \
            }   \
        }   \
    }   \
    /*DEBUG*/   \
    /*simd_abpoa_print_ag_matrix(score_t);    
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);*/   \
    simd_abpoa_ag_local_get_cigar(score_t);     \
    simd_abpoa_free_var; \
    b_score = best_score;  \
}

#define simd_abpoa_lg_var(score_t, sp, SIMDSetOne) \
    abpoa_graph_t *graph = ab->abg; abpoa_simd_matrix_t *abm = ab->abm; \
    if (abpt->gap_open != 0) err_fatal_simple("Gap open != 0 !\n");	\
    int matrix_row_n = graph->node_n, matrix_col_n = qlen + 1;	\
    int **pre_index, *pre_n, pre_i;	\
    int i, j, k, *dp_beg, *dp_beg_sn, *dp_end, *dp_end_sn, node_id, index_i, q_i;	\
    int beg, end, beg_sn, end_sn, _beg_sn, _end_sn, pre_beg_sn, pre_beg, pre_end, pre_end_sn, sn_i; \
    int pn , qp_sn, dp_sn, size; /* pn: value per SIMDi, qp_sn/dp_sn/d_sn: segmented length	*/  \
    SIMDi *DP_H, *dp_h, *pre_dp_h, *dp_f, *qp, *qi; \
    score_t *_dp_h, *_dp_f, *_qi, best_score = sp.inf_min, inf_min = sp.inf_min;	\
    int *mat = abpt->mat, best_i = 0, best_j = 0, max, max_i; score_t gap_ext = abpt->gap_ext;	\
    int w = abpt->bw < 0 ? qlen : abpt->bw;	\
    SIMDi zero = SIMDSetZeroi(), min_inf = SIMDSetOne(inf_min), gap_e = SIMDSetOne(abpt->gap_ext);	\
{   /* allocate memory */	\
    pn = sp.num_of_value, size = sp.size; 	\
    qp_sn = dp_sn = (matrix_col_n + pn - 1) / pn;	\
    qp = abm->s_mem;    \
    for (i = 0; i < qp_sn * abpt->m; ++i) qp[i] = min_inf; /* generate the query profile */	\
    for (k = 0; k < abpt->m; ++k) { /* SIMD parallelization */	\
        int *p = &mat[k * abpt->m];	\
        score_t *_qp = (score_t*)(qp + k * qp_sn);	\
        _qp[0] = 0; 	\
        for (j = 0; j < qlen; ++j) _qp[j+1] = (score_t)p[query[j]];	\
        for (j = qlen+1; j < qp_sn * pn; ++j) _qp[j] = 0;	\
    }	\
    DP_H = qp + qp_sn * abpt->m; /* H|E|H|E	*/  \
    dp_f = DP_H + dp_sn * matrix_row_n;	\
    if (abpt->use_ada || abpt->align_mode == ABPOA_EXTEND_MODE || abpt->align_mode == ABPOA_LOCAL_MODE) { /* query index */   \
        qi = dp_f + dp_sn; _qi = (score_t*)qi;    \
        for (i = 0; i <= qlen; ++i) _qi[i] = i; \
        for (i = qlen+1; i < (qlen/pn+1) * pn; ++i) _qi[i] = -1;    \
    }   \
    dp_beg = abm->dp_beg, dp_end = abm->dp_end, dp_beg_sn = abm->dp_beg_sn , dp_end_sn = abm->dp_end_sn;    \
    if (abpt->bw >= 0) {/* init min/max_rank */	\
        for (i = 1; i < graph->node_n; ++i) {	\
            node_id = graph->index_to_node_id[i];	\
            graph->node_id_to_min_rank[node_id] = graph->node_n;	\
            graph->node_id_to_max_rank[node_id] = 0;	\
        }	\
    }   \
    /* index of pre-node */	\
    pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));	\
    pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));	\
    for (i = 0; i < graph->node_n; ++i) {	\
        node_id = abpoa_graph_index_to_node_id(graph, i); /* i: node index */	\
        pre_n[i] = graph->node[node_id].in_edge_n;	\
        pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));	\
        for (j = 0; j < pre_n[i]; ++j) {	\
            pre_index[i][j] = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);	\
        }	\
    }	\
}

#define simd_abpoa_lg_first_row(score_t) { \
    /* fill the first row */	\
    if (abpt->bw >= 0) {    \
        graph->node_id_to_min_rank[0] = graph->node_id_to_max_rank[0] = 0;	\
        for (i = 0; i < graph->node[0].out_edge_n; ++i) { /* set min/max rank for next_id */	\
            int out_id = graph->node[0].out_id[i];	\
            graph->node_id_to_min_rank[out_id] = graph->node_id_to_max_rank[out_id] = 1;	\
        }   \
        dp_beg[0] = GET_DP_BEGIN(graph, w, 0, qlen), dp_end[0] = GET_DP_END(graph, w, 0, qlen);	\
    }else {    \
        dp_beg[0] = 0, dp_end[0] = qlen;	\
    }   \
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;	\
    dp_h = DP_H; _dp_h = (score_t*)dp_h;	\
    _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);	\
}

#define simd_abpoa_lg_first_dp(score_t) {   \
    simd_abpoa_lg_first_row(score_t);   \
    for (i = 0; i <= _end_sn; ++i) {	\
        dp_h[i] = min_inf;	\
    }	\
    for (i = 0; i <= dp_end[0]; ++i) { /* no SIMD parallelization */	\
        _dp_h[i] = -gap_ext * i;	\
    }	\
}

#define simd_abpoa_lg_local_first_dp(score_t) { \
    simd_abpoa_lg_first_row(score_t);   \
    for (i = 0; i <= _end_sn; ++i) { /* SIMD parallelization */	\
        dp_h[i] = zero;	\
    }\
}

#define simd_abpoa_lg_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub) {     \
    node_id = abpoa_graph_index_to_node_id(graph, index_i);	\
    SIMDi *q = qp + graph->node[node_id].base * qp_sn;	\
    dp_h = &DP_H[index_i * dp_sn]; _dp_h = (score_t*)dp_h, _dp_f = (score_t*)dp_f;	\
	\
    if (abpt->bw < 0) beg = dp_beg[index_i] = 0, end = dp_end[index_i] = qlen;	\
    else beg = dp_beg[index_i] = GET_DP_BEGIN(graph, w, node_id, qlen), end = dp_end[index_i] = GET_DP_END(graph, w, node_id, qlen);	\
    beg_sn = dp_beg_sn[index_i] = (dp_beg[index_i])/pn; end_sn = dp_end_sn[index_i] = (dp_end[index_i])/pn;	\
	\
    /* loop query */	\
    /* init h, e */	\
    _beg_sn = MAX_OF_TWO(beg_sn-1, 0); _end_sn = MIN_OF_TWO(end_sn+1, dp_sn-1);	\
    for (i = _beg_sn; i <= _end_sn; ++i) { /* SIMD parallelization */	\
        dp_h[i] = min_inf; dp_f[i] = min_inf;	\
    }	\
	\
    /* get max m and e */	\
    for (i = 0; i < pre_n[index_i]; ++i) {	\
        pre_i = pre_index[index_i][i];	\
        pre_dp_h = DP_H + pre_i * dp_sn;	\
        pre_beg = dp_beg[pre_i]; pre_end = dp_end[pre_i];	\
        pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];	\
	\
        /* set M from (pre_i, q_i-1), E from (pre_i, q_i) */	\
        _beg_sn = MAX_OF_THREE(0, pre_beg_sn, (beg-1)/pn); _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);	\
        SIMDi first = SIMDShiftRight(min_inf, SIMDTotalBytes-SIMDShiftOneN);	\
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */	\
            SIMDi remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);	\
            dp_h[sn_i] = SIMDMax(SIMDAdd(SIMDOri(first, remain), q[sn_i]), SIMDMax(SIMDSub(pre_dp_h[sn_i], gap_e), dp_h[sn_i]));	\
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);	\
        } /* now we have max(h,e) stored at dp_h */	\
    }	\
	\
    for (q_i = beg+1; q_i <= end; ++q_i) {	\
        _dp_f[q_i] = MAX_OF_TWO(_dp_h[q_i-1], _dp_f[q_i-1]) - gap_ext;	\
    }	\
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */	\
        dp_h[sn_i] = SIMDMax(dp_f[sn_i], dp_h[sn_i]);	\
    }	\
}

#define simd_abpoa_lg_local_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub) { \
    node_id = abpoa_graph_index_to_node_id(graph, index_i);	\
    SIMDi *q = qp + graph->node[node_id].base * qp_sn;	\
    dp_h = DP_H + index_i * dp_sn;	\
    _dp_h = (score_t*)dp_h, _dp_f = (score_t*)dp_f;	\
	\
    if (abpt->bw < 0) beg = dp_beg[index_i] = 0, end = dp_end[index_i] = qlen;	\
    else beg = dp_beg[index_i] = GET_DP_BEGIN(graph, w, node_id, qlen), end = dp_end[index_i] = GET_DP_END(graph, w, node_id, qlen);	\
    beg_sn = dp_beg_sn[index_i] = (dp_beg[index_i])/pn; end_sn = dp_end_sn[index_i] = (dp_end[index_i])/pn;	\
	\
    /* loop query */	\
    /* init h, e */	\
    _beg_sn = MAX_OF_TWO(beg_sn-1, 0); _end_sn = MIN_OF_TWO(end_sn+1, dp_sn-1);	\
    for (i = _beg_sn; i <= _end_sn; ++i) { /* SIMD parallelization */	\
        dp_h[i] =  zero; dp_f[i] = zero;	\
    }	\
	\
    /* get max m and e */	\
    for (i = 0; i < pre_n[index_i]; ++i) {	\
        pre_i = pre_index[index_i][i];	\
        pre_dp_h = DP_H + pre_i * dp_sn;	\
        pre_beg = dp_beg[pre_i]; pre_end = dp_end[pre_i];	\
        pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];	\
	\
        /* set M from (pre_i, q_i-1), E from (pre_i, q_i) */	\
        _beg_sn = MAX_OF_THREE(0, pre_beg_sn, (beg-1)/pn); _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);	\
        SIMDi first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);	\
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */	\
            SIMDi remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);	\
            dp_h[sn_i] = SIMDMax(SIMDAdd(SIMDOri(first, remain), q[sn_i]), SIMDMax(SIMDSub(pre_dp_h[sn_i], gap_e), dp_h[sn_i]));	\
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);	\
        } /* now we have max(h,e) stored at dp_h */	\
    }	\
	\
    for (q_i = beg+1; q_i <= end; ++q_i) {	\
        _dp_f[q_i] = MAX_OF_TWO(_dp_h[q_i-1], _dp_f[q_i-1]) - gap_ext;	\
    }	\
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */	\
        dp_h[sn_i] = SIMDMax(zero, SIMDMax(dp_f[sn_i], dp_h[sn_i]));	\
    }	\
}

#define simd_abpoa_lg_ada_rank(score_t) { \
    /* determine max pre_i, then determin current rank */	\
    int s, max_pre_i=-1, max_pre_id;	\
    score_t *_pre_dp_h;	\
    s = graph->node[node_id].base == query[max_i-1] ? abpt->match : -abpt->mismatch;	\
    for (k = 0; k < pre_n[index_i]; ++k) {	\
        pre_i = pre_index[index_i][k];	\
        dp_h = DP_H + index_i * dp_sn; pre_dp_h = DP_H + pre_i * dp_sn;	\
        _dp_h = (score_t*)dp_h; _pre_dp_h = (score_t*)pre_dp_h;	\
        if (_pre_dp_h[max_i-1] + s == _dp_h[max_i] || _pre_dp_h[max_i] - gap_ext == _dp_h[max_i]) {	\
            max_pre_i = pre_i;	\
            break;	\
        } 	\
    }	\
    max_pre_id = abpoa_graph_index_to_node_id(graph, max_pre_i);	\
    simd_abpoa_ada_rank_from_maxpre;    \
}

#define simd_abpoa_lg_global_get_max(score_t) {	\
    int in_id, in_index;	\
    for (i = 0; i < graph->node[ABPOA_SINK_NODE_ID].in_edge_n; ++i) {	\
        in_id = graph->node[ABPOA_SINK_NODE_ID].in_id[i];	\
        in_index = abpoa_graph_node_id_to_index(graph, in_id);	\
        dp_h = DP_H + in_index * dp_sn;	\
        _dp_h = (score_t*)dp_h;	\
        _set_max_score(best_score, best_i, best_j, _dp_h[qlen], in_index, qlen);	\
    }	\
}

#define simd_abpoa_lg_get_cigar(score_t) {	\
    if (abpt->ret_cigar) {	\
        simd_abpoa_lg_backtrack(score_t, DP_H, dp_sn, (score_t)abpt->match, (score_t)abpt->mismatch, (score_t)abpt->gap_ext, pre_index, pre_n, 0, 0, best_i, best_j, qlen, graph, query, n_cigar, graph_cigar);	\
    }	\
}

#define simd_abpoa_lg_local_get_cigar(score_t) {	\
    if (abpt->ret_cigar) {   \
        simd_abpoa_lg_local_backtrack(score_t, DP_H, dp_sn, (score_t)abpt->match, (score_t)abpt->mismatch, (score_t)abpt->gap_ext, pre_index, pre_n, 0, 0, best_i, best_j, qlen, graph, query, n_cigar, graph_cigar);   \
    }   \
}

// linear gap penalty: gap_open == 0
#define simd_abpoa_lg_global_align_sequence_with_graph_core(score_t, ab, query, qlen, abpt, n_cigar, graph_cigar, sp,   \
        SIMDSetOne, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, b_score) { \
    simd_abpoa_lg_var(score_t, sp, SIMDSetOne); \
    simd_abpoa_lg_first_dp(score_t);    \
    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {    \
        simd_abpoa_lg_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub);  \
        if (abpt->bw >= 0) {    \
            if (abpt->use_ada) { /* set current rank, set min/max rank for next nodes */    \
                simd_abpoa_max(score_t, SIMDSetIfGreater, SIMDGetIfGreater);    \
                /*DEBUG*/   \
                /*printf("adaMax: %d, %d\n", max_i, max);*/ \
                simd_abpoa_lg_ada_rank(score_t);    \
            } else { /* set min/max_rank for next nodes */ \
                simd_abpoa_rank;    \
            }   \
        }   \
    }   \
    simd_abpoa_lg_global_get_max(score_t);  \
    /*DEBUG*/   \
    /*simd_abpoa_print_lg_matrix(score_t);
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);*/   \
    simd_abpoa_lg_get_cigar(score_t);   \
    simd_abpoa_free_var; \
    b_score = best_score;  \
} 

#define simd_abpoa_lg_extend_align_sequence_with_graph_core(score_t, graph, query, qlen, abpt, n_cigar, graph_cigar, sp,    \
    SIMDSetOne, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, b_score) { \
    simd_abpoa_lg_var(score_t, sp, SIMDSetOne); \
    simd_abpoa_lg_first_dp(score_t);    \
    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {    \
        simd_abpoa_lg_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub);  \
        simd_abpoa_max(score_t, SIMDSetIfGreater, SIMDGetIfGreater);    \
        /*DEBUG*/   \
        /* printf("extendMax: %d, %d\n", max_i, max); */    \
        if (abpt->zdrop > 0 && best_score - max > abpt->zdrop + gap_ext * abs((best_i-index_i)-(best_j-max_i))) \
            break;  \
        _set_max_score(best_score, best_i, best_j, max, index_i, max_i);    \
        if (abpt->bw >= 0) {    \
            if (abpt->use_ada) { /* set current rank, set min/max rank for next nodes */  \
                simd_abpoa_lg_ada_rank(score_t);    \
            } else { /* set min/max_rank for next nodes */ \
                simd_abpoa_rank; \
            }   \
        }   \
    }   \
    /*DEBUG*/   \
    /* simd_abpoa_print_lg_matrix(score_t);
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score); */  \
    simd_abpoa_lg_get_cigar(score_t);   \
    simd_abpoa_free_var; \
    b_score = best_score;  \
}

#define simd_abpoa_lg_local_align_sequence_with_graph_core(score_t, graph, query, qlen, abpt, n_cigar, graph_cigar, sp,     \
    SIMDSetOne, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, b_score) { \
    simd_abpoa_lg_var(score_t, sp, SIMDSetOne); \
    simd_abpoa_lg_local_first_dp(score_t);  \
    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {    \
        simd_abpoa_lg_local_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub);    \
        simd_abpoa_max(score_t, SIMDSetIfGreater, SIMDGetIfGreater);    \
        /*DEBUG*/   \
        /* printf("localMax: %d, %d\n", max_i, max); */     \
        _set_max_score(best_score, best_i, best_j, max, index_i, max_i);    \
        if (abpt->bw >= 0) { \
            if (max > 0 && abpt->use_ada) { /* set current rank, set min/max rank for next nodes */    \
                simd_abpoa_lg_ada_rank(score_t);    \
            } else { /* set min/max_rank for next nodes */  \
                simd_abpoa_rank;        \
            }   \
        }   \
    }   \
    /*DEBUG*/   \
    /* simd_abpoa_print_lg_matrix(score_t);
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score); */ \
    simd_abpoa_lg_local_get_cigar(score_t); \
    simd_abpoa_free_var; \
    b_score = best_score;  \
} 

abpoa_simd_matrix_t *abpoa_init_simd_matrix(void) {
    abpoa_simd_matrix_t *abm = (abpoa_simd_matrix_t*)_err_malloc(sizeof(abpoa_simd_matrix_t));
    abm->s_msize = 0; abm->s_mem = NULL;
    abm->rang_m = 0; abm->dp_beg = NULL; abm->dp_end = NULL; abm->dp_beg_sn = NULL; abm->dp_end_sn = NULL;
    return abm;
}

void abpoa_free_simd_matrix(abpoa_simd_matrix_t *abm) {
    if (abm->s_mem) free(abm->s_mem);
    if (abm->dp_beg) {
        free(abm->dp_beg);
        free(abm->dp_end);
        free(abm->dp_beg_sn);
        free(abm->dp_end_sn);
    }
    free(abm);
}

// realloc memory everytime the graph is updated (nodes are updated already)
// * index_to_node_id/node_id_to_index/node_id_to_max/min_rank/remain
// * qp, DP_HE/H (if ag/lg), dp_f, qi (if ada/extend/local), hd (if ada), backtrack_z (if cigar)
// * dp_beg/end, dp_beg/end_sn if band
// * pre_n, pre_index
void simd_abpoa_realloc(abpoa_t *ab, int qlen, abpoa_para_t *abpt, SIMD_para_t sp) {
    int pn = sp.num_of_value, size = sp.size, sn = (qlen + sp.num_of_value) / pn, node_n = ab->abg->node_n;
    int s_msize = sn * abpt->m * size; // qp
    if (abpt->gap_open == 0) s_msize += (sn * (node_n + 1) * size); // DP_H + dp_f
    else s_msize += (sn * (2 * node_n + 1) * size); // DP_HE + dp_f
    if (abpt->use_ada || abpt->align_mode == ABPOA_EXTEND_MODE || abpt->align_mode == ABPOA_LOCAL_MODE) // qi
        s_msize += sn * size;
    if (abpt->gap_open != 0) {
        if (abpt->ret_cigar) // hd + backtrack_z
            s_msize += sn * (node_n - 1) * size;
        else if (abpt->use_ada) // hd
            s_msize += sn * size;
    }

    if (s_msize > ab->abm->s_msize) {
        if (ab->abm->s_mem) free(ab->abm->s_mem);
        ab->abm->s_msize = MAX_OF_TWO(ab->abm->s_msize << 1, s_msize);
        ab->abm->s_mem = (SIMDi*)SIMDMalloc(ab->abm->s_msize, size);
    }

    if (node_n > ab->abm->rang_m) {
        ab->abm->rang_m = MAX_OF_TWO(ab->abm->rang_m << 1, node_n);
        ab->abm->dp_beg = (int*)_err_realloc(ab->abm->dp_beg, ab->abm->rang_m * sizeof(int));
        ab->abm->dp_end = (int*)_err_realloc(ab->abm->dp_end, ab->abm->rang_m * sizeof(int));
        ab->abm->dp_beg_sn = (int*)_err_realloc(ab->abm->dp_beg_sn, ab->abm->rang_m * sizeof(int));
        ab->abm->dp_end_sn = (int*)_err_realloc(ab->abm->dp_end_sn, ab->abm->rang_m * sizeof(int));
    }
}

// based on min,max and w
// set score_n and id_n (8, 16, 32 bits)
int simd_abpoa_align_sequence_with_graph(abpoa_t *ab, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar) {
    if (abpt->simd_flag == 0) return ada_abpoa_banded_global_align_sequence_with_graph(ab->abg, query, qlen, abpt, n_cigar, graph_cigar);

    int max_score;
    if (abpt->simd_flag & SIMD_AVX512F && !(abpt->simd_flag & SIMD_AVX512BW)) max_score = 32768; // AVX512F has no 8/16 bits operations
    else {
        int node_n = ab->abg->node_n;
        max_score = MAX_OF_THREE(qlen * abpt->match, MIN_OF_TWO(qlen, node_n) * abpt->mismatch + abpt->gap_open + abpt->gap_ext * abs(qlen-node_n), abpt->gap_open + node_n * abpt->gap_ext + abpt->gap_open + qlen * abpt->gap_ext);
    }
    int bits = 0;
    if (max_score <= 127 - abpt->mismatch - abpt->gap_open - abpt->gap_ext) { // DP_H/E/F: 8  bits
        _simd_p8.inf_min = MAX_OF_TWO(INT8_MIN + abpt->mismatch, INT8_MIN + abpt->gap_open + abpt->gap_ext);
        simd_abpoa_realloc(ab, qlen, abpt, _simd_p8);
        bits = 8;
    } else if (max_score <= 32767 - abpt->mismatch - abpt->gap_open - abpt->gap_ext) { // DP_H/E/F: 16 bits
        _simd_p16.inf_min = MAX_OF_TWO(INT16_MIN + abpt->mismatch, INT16_MIN + abpt->gap_open + abpt->gap_ext);
        simd_abpoa_realloc(ab, qlen, abpt, _simd_p16);
        bits = 16;
    } else { 
        _simd_p32.inf_min = MAX_OF_TWO(INT32_MIN + abpt->mismatch, INT32_MIN + abpt->gap_open + abpt->gap_ext);
        simd_abpoa_realloc(ab, qlen, abpt, _simd_p32);
        bits = 32;
    }


    int _best_score=0;
    if (abpt->align_mode == ABPOA_GLOBAL_MODE) {
        if (bits == 8) {
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_global_align_sequence_with_graph_core(int8_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p8, \
                        SIMDSetOnei8, SIMDMaxi8, SIMDAddi8, SIMDSubi8, SIMDShiftOneNi8, SIMDSetIfGreateri8, SIMDGetIfGreateri8, _best_score);
            } else {
                simd_abpoa_lg_global_align_sequence_with_graph_core(int8_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p8, \
                        SIMDSetOnei8, SIMDMaxi8, SIMDAddi8, SIMDSubi8, SIMDShiftOneNi8, SIMDSetIfGreateri8, SIMDGetIfGreateri8,  _best_score);
            }
        } else if (bits == 16) {
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_global_align_sequence_with_graph_core(int16_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p16, \
                        SIMDSetOnei16, SIMDMaxi16, SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, _best_score);
            } else {
                simd_abpoa_lg_global_align_sequence_with_graph_core(int16_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p16, \
                        SIMDSetOnei16, SIMDMaxi16, SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, _best_score);
            }
        } else { // 2147483647, DP_H/E/F: 32 bits
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_global_align_sequence_with_graph_core(int32_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p32, \
                        SIMDSetOnei32, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, _best_score);
            } else {
                simd_abpoa_lg_global_align_sequence_with_graph_core(int32_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p32,    \
                        SIMDSetOnei32, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, _best_score);
            }
        }
    } else if (abpt->align_mode == ABPOA_LOCAL_MODE) {
        if (bits == 8) {
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_local_align_sequence_with_graph_core(int8_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p8, \
                        SIMDSetOnei8, SIMDMaxi8, SIMDAddi8, SIMDSubi8, SIMDShiftOneNi8, SIMDSetIfGreateri8, SIMDGetIfGreateri8, _best_score);
            } else {
                simd_abpoa_lg_local_align_sequence_with_graph_core(int8_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p8, \
                        SIMDSetOnei8, SIMDMaxi8, SIMDAddi8, SIMDSubi8, SIMDShiftOneNi8, SIMDSetIfGreateri8, SIMDGetIfGreateri8, _best_score);
            }
        } else if (bits == 16) {
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_local_align_sequence_with_graph_core(int16_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p16, \
                        SIMDSetOnei16, SIMDMaxi16, SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, _best_score);
            } else {
                simd_abpoa_lg_local_align_sequence_with_graph_core(int16_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p16, \
                        SIMDSetOnei16, SIMDMaxi16, SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, _best_score);
            }
        } else { // 2147483647, DP_H/E/F: 32 bits
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_local_align_sequence_with_graph_core(int32_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p32, \
                        SIMDSetOnei32, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, _best_score);
            } else {
                simd_abpoa_lg_local_align_sequence_with_graph_core(int32_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p32,    \
                        SIMDSetOnei32, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32,  _best_score);
            }
        }
    } else if (abpt->align_mode == ABPOA_EXTEND_MODE) {
        if (bits == 8) {
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_extend_align_sequence_with_graph_core(int8_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p8, \
                        SIMDSetOnei8, SIMDMaxi8, SIMDAddi8, SIMDSubi8, SIMDShiftOneNi8, SIMDSetIfGreateri8, SIMDGetIfGreateri8, _best_score);
            } else {
                simd_abpoa_lg_extend_align_sequence_with_graph_core(int8_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p8, \
                        SIMDSetOnei8, SIMDMaxi8, SIMDAddi8, SIMDSubi8, SIMDShiftOneNi8, SIMDSetIfGreateri8, SIMDGetIfGreateri8, _best_score);
            }
        } else if (bits == 16) {
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_extend_align_sequence_with_graph_core(int16_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p16, \
                        SIMDSetOnei16, SIMDMaxi16, SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, _best_score);
            } else {
                simd_abpoa_lg_extend_align_sequence_with_graph_core(int16_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p16, \
                        SIMDSetOnei16, SIMDMaxi16, SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, _best_score);
            }
        } else { // 2147483647, DP_H/E/F: 32 bits
            if (abpt->gap_open > 0) {
                simd_abpoa_ag_extend_align_sequence_with_graph_core(int32_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p32, \
                        SIMDSetOnei32, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, _best_score);
            } else {
                simd_abpoa_lg_extend_align_sequence_with_graph_core(int32_t, ab, query, qlen, abpt, n_cigar, graph_cigar, _simd_p32,    \
                        SIMDSetOnei32, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, _best_score);
            }
        }
    } else if (abpt->align_mode == ABPOA_SEMI_MODE) { // TODO semi-global
    } else err_fatal_core(__func__, "Unknown align mode. (%d)\n", abpt->align_mode);
    return _best_score;
}
