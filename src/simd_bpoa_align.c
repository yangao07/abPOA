#include <stdio.h>
#include <stdlib.h>
#include "bpoa_graph.h"
#include "bpoa_align.h"
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

#define print_simd(s, str, score_t) {                           \
    int i; score_t *_a = (score_t*)s;                           \
    printf("%s\t", str);                                        \
    for (i = 0; i < SIMDTotalBytes / (int)sizeof(score_t); ++i) {                                  \
        printf("%d\t", _a[i]);                                  \
    } printf("\n");                                             \
}

#define simd_bpoa_backtrack(score_t, DP_HE, dp_sn, match, mis, gap_e, pre_index, pre_n, backtrack_z, d_sn, best_i, best_j, graph, query, n_cigar, graph_cigar) { \
    int i, j, k, pre_i;                                                                                                     \
    SIMDi *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *bz;                                                                          \
    score_t *_dp_h, *_pre_dp_h, *_dp_e, *_pre_dp_e;                                                                         \
    int n_c = 0, s, m_c = 0, id, which, last_which;                                                                         \
    int op_shift[3] = {HOP_OFF_SET, EOP_OFF_SET, FOP_OFF_SET};                                                              \
    score_t d, *_bz; bpoa_cigar_t *cigar = 0;                                                                               \
                                                                                                                            \
    i = best_i, j = best_j, id = bpoa_graph_index_to_node_id(graph, i), last_which = 0;                                     \
    while (i > 0 && j > 0) {                                                                                                \
        bz = backtrack_z + (i-1) * d_sn;                                                                                    \
        _bz = (score_t*)bz;                                                                                                 \
        d = _bz[j];                                                                                                         \
        which = (d >> op_shift[last_which]) & 3;                                                                            \
        if (which == 0) { /* match */                                                                                       \
            cigar = bpoa_push_cigar(&n_c, &m_c, cigar, POA_CMATCH, 1, id, j-1);                                             \
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
            id = bpoa_graph_index_to_node_id(graph, i);                                                                     \
            j--; last_which = which;                                                                                        \
        } else if (which == 1) { /* deletion */                                                                             \
            cigar = bpoa_push_cigar(&n_c, &m_c, cigar, POA_CDEL, 1, id, j-1);                                               \
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
            id = bpoa_graph_index_to_node_id(graph, i);                                                                     \
            last_which = which;                                                                                             \
        } else { /* insertion */                                                                                            \
            cigar = bpoa_push_cigar(&n_c, &m_c, cigar, POA_CINS, 1, id, j-1);                                               \
            j--; last_which = which;                                                                                        \
        }                                                                                                                   \
    }                                                                                                                       \
    if (j > 0) cigar = bpoa_push_cigar(&n_c, &m_c, cigar, POA_CSOFT_CLIP, j, -1, j-1);                                      \
    /* reverse cigar */                                                                                                     \
    *graph_cigar = bpoa_reverse_cigar(n_c, cigar);                                                                          \
    *n_cigar = n_c;                                                                                                         \
    /* DEBUG */                                                                                                             \
    bpoa_print_cigar(n_c, *graph_cigar, graph);                                                                             \
}

#define simd_bpoa_banded_global_align_sequence_with_graph_core(graph, query, qlen, bpt, n_cigar, graph_cigar, sp, score_t,  \
    SIMDSetOne, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, SIMDSetIfEqual, b_score) {    \
    int **pre_index, *pre_n, pre_i;                                                                                         \
    int target_node_n = graph->node_n - 2, matrix_row_n = graph->node_n, matrix_col_n = qlen + 1;                           \
    int i, j, k, w, *dp_beg, *dp_beg_sn, *dp_end, *dp_end_sn, node_id, index_i, q_i;                                        \
    int pn , qp_sn,d_sn, dp_sn, size; /* pn: value per SIMDi, qp_sn/dp_sn/d_sn: segmented length */                         \
                                                                                                                            \
    SIMDi *qp, gap_e, gap_oe, min_inf;                                                                                      \
    SIMDi *DP_HE, *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *dp_f, tmp_e;                                                         \
    SIMDi *backtrack_z, *z; /* backtrack cell: f<<4|e<<2|h, MATCH:0, DELETION:1, INSERTION:2 */                             \
    SIMDi *hd, fd, ed, hm, he, hf, ee, em, ff, fm; uint8_t m0=0x0, e1=0x1, f2=0x2;                                          \
                                                                                                                            \
    score_t *_dp_h, *_dp_e, *_dp_f, best_score = sp.inf_min, inf_min = sp.inf_min; int best_i=0, best_j=0;                  \
    int *mat = bpt->mat, gap_open = bpt->gap_open, gap_ext = bpt->gap_ext, gap_OE = bpt->gap_open + bpt->gap_ext;           \
                                                                                                                            \
    {   /* allocate memory */                                                                                               \
        pn = sp.num_of_value, size = sp.size;                                                                               \
        qp_sn = d_sn = dp_sn = (matrix_col_n + pn - 1) / pn;                                                                \
        qp = (SIMDi*)SIMDMalloc((qp_sn * bpt->m + dp_sn * (2 * matrix_row_n + 1)) * size, size);                            \
        DP_HE = qp + qp_sn * bpt->m; /* H|E|H|E */                                                                          \
        dp_f = DP_HE + 2 * dp_sn * matrix_row_n;                                                                            \
                                                                                                                            \
        dp_beg = (int*)_err_malloc((matrix_row_n) * sizeof(int)); dp_end = (int*)_err_malloc(matrix_row_n * sizeof(int));   \
        dp_beg_sn = (int*)_err_malloc((matrix_row_n) * sizeof(int)); dp_end_sn = (int*)_err_malloc(matrix_row_n * sizeof(int)); \
        /* when w <= 0, do whole global */                                                                                  \
        w = bpt->bw <= 0 ? qlen : bpt->bw;                                                                                  \
                                                                                                                            \
        min_inf = SIMDSetOne(inf_min); gap_e = SIMDSetOne(bpt->gap_ext), gap_oe = SIMDSetOne(bpt->gap_open + bpt->gap_ext); \
        if (graph_cigar && n_cigar) {                                                                                       \
            /* for backtrack */                                                                                             \
            hd = (SIMDi*)SIMDMalloc(d_sn * (target_node_n + 1) * size, size);                                               \
            backtrack_z = hd + d_sn;                                                                                        \
                                                                                                                            \
            hm = SIMDSetOne(m0 << HOP_OFF_SET), he = SIMDSetOne(e1 << HOP_OFF_SET), hf = SIMDSetOne(f2 << HOP_OFF_SET);     \
            em = SIMDSetOne(m0 << EOP_OFF_SET), ee = SIMDSetOne(e1 << EOP_OFF_SET), fm = SIMDSetOne(m0 << FOP_OFF_SET), ff = SIMDSetOne(f2 << FOP_OFF_SET); \
        }                                                                                                                   \
                                                                                                                            \
        pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));                                                       \
        pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));                                                            \
    }                                                                                                                       \
                                                                                                                            \
    {   /* generate the query profile */                                                                                    \
        for (k = 0; k < bpt->m; ++k) { /* SIMD parallelization */                                                           \
            int *p = &mat[k * bpt->m];                                                                                      \
            score_t *_qp = (score_t*)(qp + k * qp_sn);                                                                      \
            _qp[0] = 0;                                                                                                     \
            for (j = 0; j < qlen; ++j) _qp[j+1] = (score_t)p[query[j]];                                                     \
        }                                                                                                                   \
                                                                                                                            \
        /* index of pre-node */                                                                                             \
        for (i = 0; i < graph->node_n; ++i) {                                                                               \
            node_id = bpoa_graph_index_to_node_id(graph, i); /* i: node index */                                            \
            pre_n[i] = graph->node[node_id].in_edge_n;                                                                      \
            pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));                                                       \
            for (j = 0; j < pre_n[i]; ++j) {                                                                                \
                pre_index[i][j] = bpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);                        \
            }                                                                                                               \
        }                                                                                                                   \
    }                                                                                                                       \
                                                                                                                            \
    {   /* DP loop */                                                                                                       \
        /* DP cell: H[i,j], E[i+1,j], F[i,j+1] */                                                                           \
        int beg, end, beg_sn, end_sn, _end_sn, pre_beg_sn, pre_end, pre_end_sn, sn_i;                                       \
                                                                                                                            \
        /* fill the first row */                                                                                            \
        dp_beg[0] = GET_DP_BEGIN(graph, w, 0); dp_end[0] = GET_DP_END(graph, w, 0);                                         \
        dp_beg_sn[0] = (dp_beg[0]+pn)/pn-1; dp_end_sn[0] = (dp_end[0]+pn)/pn-1;                                             \
                                                                                                                            \
        dp_h = DP_HE, dp_e = dp_h + dp_sn;                                                                                  \
        _dp_h = (score_t*)dp_h, _dp_e = (score_t*)dp_e;                                                                     \
        _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);                                                                      \
        for (i = 0; i <= _end_sn; ++i) {                                                                                    \
            SIMDStorei(&dp_h[i], min_inf); SIMDStorei(&dp_e[i], min_inf);                                                   \
        }                                                                                                                   \
        for (i = 1; i <= dp_end[0]; ++i) { /* no SIMD parallelization */                                                    \
            _dp_h[i] = -gap_open - gap_ext * i;                                                                             \
        } _dp_h[0] = 0; _dp_e[0] = -(gap_OE);                                                                               \
                                                                                                                            \
                                                                                                                            \
        for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {                                                            \
            node_id = bpoa_graph_index_to_node_id(graph, index_i);                                                          \
                                                                                                                            \
            SIMDi *q = qp + graph->node[node_id].base * qp_sn, a, b, c, d;                                                  \
            dp_h = DP_HE + index_i * 2 * dp_sn, dp_e = dp_h + dp_sn;                                                        \
            _dp_h = (score_t*)dp_h, _dp_e = (score_t*)dp_e, _dp_f = (score_t*)dp_f;                                         \
                                                                                                                            \
            beg = dp_beg[index_i] = GET_DP_BEGIN(graph, w, index_i); end = dp_end[index_i] = GET_DP_END(graph, w, index_i); \
            beg_sn = dp_beg_sn[index_i] = (dp_beg[index_i]+pn)/pn-1; end_sn = dp_end_sn[index_i] = (dp_end[index_i]+pn)/pn-1; \
                                                                                                                            \
            /* loop query */                                                                                                \
            /* initialize h, e */                                                                                           \
            _end_sn = MIN_OF_TWO(end_sn+1, dp_sn-1);                                                                        \
            for (i = beg_sn; i <= end_sn; ++i) { /* SIMD parallelization */                                                 \
                SIMDStorei(&dp_h[i], min_inf); SIMDStorei(&dp_e[i], min_inf); SIMDStorei(&dp_f[i], min_inf);                \
            }                                                                                                               \
                                                                                                                            \
            for (i = 0; i < pre_n[index_i]; ++i) {                                                                          \
                pre_i = pre_index[index_i][i];                                                                              \
                pre_dp_h = DP_HE + pre_i * 2 * dp_sn, pre_dp_e = pre_dp_h + dp_sn;                                          \
                pre_end = dp_end[pre_i];                                                                                    \
                pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                               \
                /* set M from (pre_i, q_i-1) */                                                                             \
                _end_sn = MIN_OF_TWO((pre_end + pn)/pn-1, dp_sn-1);                                                         \
                SIMDi first, remain;                                                                                        \
                first = SIMDShiftRight(min_inf, SIMDTotalBytes-SIMDShiftOneN);                                              \
                for (sn_i = pre_beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                               \
                    a = SIMDLoadi(&pre_dp_h[sn_i]); b = SIMDLoadi(&dp_h[sn_i]);                                             \
                    remain = SIMDShiftLeft(a, SIMDShiftOneN);                                                               \
                    SIMDStorei(&dp_h[sn_i], SIMDMax(SIMDOri(first, remain), b));                                            \
                    first = SIMDShiftRight(a, SIMDTotalBytes-SIMDShiftOneN);                                                \
                }                                                                                                           \
                /* set E from (pre_i, q_i) */                                                                             \
                for (sn_i = pre_beg_sn; sn_i <= pre_end_sn; ++sn_i) { /* SIMD parallelization */                            \
                    a = SIMDLoadi(&pre_dp_e[sn_i]); b = SIMDLoadi(&dp_e[sn_i]);                                             \
                    SIMDStorei(&dp_e[sn_i], SIMDMax(a, b));                                                                 \
                }                                                                                                           \
            }                                                                                                               \
                                                                                                                            \
            if (n_cigar && graph_cigar) {                                                                                   \
                z = backtrack_z + (index_i-1) * d_sn;                                                                       \
                /* compare M, E, and F */                                                                                   \
                for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */                                    \
                    a = SIMDLoadi(&dp_h[sn_i]); b = SIMDLoadi(&q[sn_i]); c = SIMDLoadi(&dp_e[sn_i]); d = SIMDLoadi(&hd[sn_i]); \
                    a = SIMDAdd(a, b); /* q[0] == 0 */                                                                      \
                    SIMDGetIfGreater(d, a, c, a, he, hm);                                                                   \
                    SIMDStorei(hd+sn_i, d); SIMDStorei(&dp_h[sn_i], a);                                                     \
                }                                                                                                           \
                                                                                                                            \
                /* set F from (index_i, q_i-1) */                                                                           \
                for (q_i = beg+1; q_i <= end; ++q_i) { /* no SIMD parallelization */                                        \
                    _dp_f[q_i] = MAX_OF_TWO(_dp_h[q_i-1] - gap_OE, _dp_f[q_i-1] - gap_ext);                                 \
                }                                                                                                           \
                                                                                                                            \
                for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */                                    \
                    a = SIMDLoadi(&dp_h[sn_i]); b = SIMDLoadi(&dp_f[sn_i]); c = SIMDLoadi(&dp_e[sn_i]); d = SIMDLoadi(&hd[sn_i]); \
                    /* h, hd for current cell */                                                                            \
                    SIMDGetIfGreater(d, a, b, a, hf, d);                                                                    \
                    SIMDStorei(&dp_h[sn_i], a);\
                                                                                                                            \
                    /* e, ed for next cell */                                                                               \
                    tmp_e = SIMDSub(a, gap_oe);                                                                             \
                    c = SIMDSub(c, gap_e);                                                                                  \
                    if (gap_open) {                                                                                         \
                        SIMDGetIfGreater(ed, c, c, tmp_e, ee, em);                                                          \
                        SIMDStorei(&dp_e[sn_i], c);                                                                         \
                    } else {                                                                                                \
                        SIMDSetIfEqual(ed, d, hm, em, ee);                                                                  \
                        SIMDStorei(&dp_e[sn_i], tmp_e);                                                                     \
                    }                                                                                                       \
                                                                                                                            \
                    /* fd for next cell */                                                                                  \
                    b = SIMDSub(b, gap_e);                                                                                  \
                    if (gap_open) {                                                                                         \
                        SIMDSetIfGreater(fd, b, tmp_e, ff, fm);                                                             \
                    } else {                                                                                                \
                        SIMDSetIfEqual(fd, d, hm, fm ,ff);                                                                  \
                    }                                                                                                       \
                    SIMDStorei(&z[sn_i], SIMDOri(SIMDOri(d, ed), fd));                                                      \
                }                                                                                                           \
            } else {                                                                                                        \
                /* compare M, E, and F */                                                                                   \
                for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */                                    \
                    a = SIMDLoadi(&dp_h[sn_i]); b = SIMDLoadi(&q[sn_i]); c = SIMDLoadi(&dp_e[sn_i]);                        \
                    a = SIMDAdd(a, b); /* q[0] == 0 */                                                                      \
                    SIMDStorei(&dp_h[sn_i], SIMDMax(a, c));                                                                 \
                }                                                                                                           \
                                                                                                                            \
                /* set F from (index_i, q_i-1) */                                                                           \
                _dp_f[beg] = inf_min;                                                                                       \
                for (q_i = beg+1; q_i <= end; ++q_i) { /* no SIMD parallelization */                                        \
                    _dp_f[q_i] = MAX_OF_TWO(_dp_h[q_i-1] - gap_OE, _dp_f[q_i-1] - gap_ext);                                 \
                }                                                                                                           \
                                                                                                                            \
                for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */                                    \
                    a = SIMDLoadi(&dp_h[sn_i]); b = SIMDLoadi(&dp_f[sn_i]); c = SIMDLoadi(&dp_e[sn_i]);                     \
                    /* h, hd for current cell */                                                                            \
                    a = SIMDMax(a, b);                                                                                      \
                    SIMDStorei(&dp_h[sn_i], a);                                                                             \
                                                                                                                            \
                    /* e for next cell */                                                                                   \
                    tmp_e = SIMDSub(a, gap_oe);                                                                             \
                    c = SIMDSub(c, gap_e);                                                                                  \
                    SIMDStorei(&dp_e[sn_i], SIMDMax(c, tmp_e));                                                             \
                }                                                                                                           \
            }                                                                                                               \
        }                                                                                                                   \
    }                                                                                                                       \
    /* DEBUG */                                                                                                             \
    for (j = 0; j <= target_node_n; ++j) {                                                                                  \
        printf("index: %d\t", j);                                                                                           \
        dp_h = DP_HE + j * 2 * dp_sn, dp_e = dp_h + dp_sn;                                                                  \
        _dp_h = (score_t*)dp_h, _dp_e = (score_t*)dp_e;                                                                     \
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {                                                                          \
            printf("%d:(%d,%d)\t", i, _dp_h[i], _dp_e[i]);                                                                  \
        } printf("\n");                                                                                                     \
    }                                                                                                                       \
    int in_id, in_index;                                                                                                    \
    for (i = 0; i < graph->node[POA_SINK_NODE_ID].in_edge_n; ++i) { /* find best backtrack position */                      \
        in_id = graph->node[POA_SINK_NODE_ID].in_id[i];                                                                     \
        in_index = bpoa_graph_node_id_to_index(graph, in_id);                                                               \
        dp_h = DP_HE + in_index * 2 * dp_sn;                                                                                \
        _dp_h = (score_t*)dp_h;                                                                                             \
        _set_max_score(best_score, best_i, best_j, _dp_h[qlen], in_index, qlen);                                            \
    }                                                                                                                       \
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);                                                     \
    { /* backtrack from best score */                                                                                       \
        if (n_cigar && graph_cigar) {                                                                                       \
            simd_bpoa_backtrack(score_t, DP_HE, dp_sn, (score_t)bpt->match, (score_t)bpt->mismatch, (score_t)bpt->gap_ext,  \
                    pre_index, pre_n, backtrack_z, d_sn, best_i, best_j, graph, query, n_cigar, graph_cigar);               \
        }                                                                                                                   \
    }                                                                                                                       \
                                                                                                                            \
    { /* free variables */                                                                                                  \
        SIMDFree(qp); free(dp_beg); free(dp_end); free(dp_beg_sn); free(dp_end_sn);                                         \
        if (n_cigar && graph_cigar) SIMDFree(hd);                                                                           \
        for (i = 0; i < graph->node_n; ++i) free(pre_index[i]);                                                             \
        free(pre_index); free(pre_n);                                                                                       \
    }                                                                                                                       \
    b_score = best_score;                                                                                                   \
}

//TODO realloc with exisiting mem
//TODO only store HE in band range
int simd16_bpoa_banded_global_align_sequence_with_graph_core(bpoa_graph_t *graph, uint8_t *query, int qlen, bpoa_para_t *bpt, int *n_cigar, bpoa_cigar_t **graph_cigar, SIMD_para_t sp) {
    int **pre_index, *pre_n, pre_i;
    int target_node_n = graph->node_n - 2, matrix_row_n = graph->node_n, matrix_col_n = qlen + 1;
    int i, j, k, w, *dp_beg, *dp_beg_sn, *dp_end, *dp_end_sn, node_id, index_i, q_i;
    int pn , qp_sn,d_sn, dp_sn, size; // pn: value per SIMDi, qp_sn/dp_sn/d_sn: segmented length

    SIMDi *qp, gap_e, gap_oe, min_inf;
    SIMDi *DP_HE, *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *dp_f, tmp_e;
    SIMDi *backtrack_z, *z; // backtrack cell: f<<4|e<<2|h, MATCH:0, DELETION:1, INSERTION:2
    SIMDi *hd, fd, ed, hm, he, hf, ee, em, ff, fm; uint8_t m0=0x0, e1=0x1, f2=0x2;
    
    int16_t *_dp_h, *_dp_e, *_dp_f, best_score = sp.inf_min, inf_min = sp.inf_min; int best_i=0, best_j=0;
    int *mat = bpt->mat;
    int16_t gap_open = bpt->gap_open, gap_ext = bpt->gap_ext, gap_OE = bpt->gap_open + bpt->gap_ext;

    {   // allocate memory 
        pn = sp.num_of_value, size = sp.size; 
        qp_sn = d_sn = dp_sn = (matrix_col_n + pn - 1) / pn;
        qp = (SIMDi*)SIMDMalloc((qp_sn * bpt->m + dp_sn * (2 * matrix_row_n + 1)) * size, size);
        DP_HE = qp + qp_sn * bpt->m; // H|E|H|E
        dp_f = DP_HE + 2 * dp_sn * matrix_row_n;

        dp_beg = (int*)_err_malloc((matrix_row_n) * sizeof(int)); dp_end = (int*)_err_malloc(matrix_row_n * sizeof(int));
        dp_beg_sn = (int*)_err_malloc((matrix_row_n) * sizeof(int)); dp_end_sn = (int*)_err_malloc(matrix_row_n * sizeof(int));
        // when w <= 0, do whole global
        w = bpt->bw <= 0 ? qlen : bpt->bw;

        min_inf = SIMDSetOnei16(inf_min); gap_e = SIMDSetOnei16(bpt->gap_ext), gap_oe = SIMDSetOnei16(bpt->gap_open + bpt->gap_ext);
        if (graph_cigar && n_cigar) {
            hd = (SIMDi*)SIMDMalloc(d_sn * (target_node_n + 1) * size, size);
            backtrack_z = hd + d_sn;

            // for backtrack
            hm = SIMDSetOnei16(m0 << HOP_OFF_SET), he = SIMDSetOnei16(e1 << HOP_OFF_SET), hf = SIMDSetOnei16(f2 << HOP_OFF_SET);
            em = SIMDSetOnei16(m0 << EOP_OFF_SET), ee = SIMDSetOnei16(e1 << EOP_OFF_SET), fm = SIMDSetOnei16(m0 << FOP_OFF_SET), ff = SIMDSetOnei16(f2 << FOP_OFF_SET);
        }

        pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
        pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));
    }

    {   // generate the query profile
        for (k = 0; k < bpt->m; ++k) { // SIMD parallelization
            int *p = &mat[k * bpt->m];
            int16_t *_qp = (int16_t*)(qp + k * qp_sn);
            _qp[0] = 0; 
            for (j = 0; j < qlen; ++j) _qp[j+1] = (int16_t)p[query[j]];
        }

        // index of pre-node
        for (i = 0; i < graph->node_n; ++i) {
            node_id = bpoa_graph_index_to_node_id(graph, i); // i: node index
            pre_n[i] = graph->node[node_id].in_edge_n;
            pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
            for (j = 0; j < pre_n[i]; ++j) {
                pre_index[i][j] = bpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
            }
        }
    }

    {   // DP loop
        // DP cell: H[i,j], E[i+1,j], F[i,j+1]
        int beg, end, beg_sn, end_sn, _end_sn, pre_beg_sn, pre_end, pre_end_sn, sn_i;

        // fill the first row
        dp_beg[0] = GET_DP_BEGIN(graph, w, 0); dp_end[0] = GET_DP_END(graph, w, 0);
        dp_beg_sn[0] = (dp_beg[0]+pn)/pn-1; dp_end_sn[0] = (dp_end[0]+pn)/pn-1;

        dp_h = DP_HE, dp_e = dp_h + dp_sn;
        _dp_h = (int16_t*)dp_h, _dp_e = (int16_t*)dp_e;
        _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);
        for (i = 0; i <= _end_sn; ++i) {
            SIMDStorei(&dp_h[i], min_inf); SIMDStorei(&dp_e[i], min_inf);
        }
        for (i = 1; i <= dp_end[0]; ++i) { // no SIMD parallelization
            _dp_h[i] = -gap_open - gap_ext * i;
        } _dp_h[0] = 0; _dp_e[0] = -(gap_OE);

        for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {
            node_id = bpoa_graph_index_to_node_id(graph, index_i);

            SIMDi *q = qp + graph->node[node_id].base * qp_sn, a, b, c, d;
            dp_h = DP_HE + index_i * 2 * dp_sn, dp_e = dp_h + dp_sn;
            _dp_h = (int16_t*)dp_h, _dp_e = (int16_t*)dp_e, _dp_f = (int16_t*)dp_f;

            beg = dp_beg[index_i] = GET_DP_BEGIN(graph, w, index_i); end = dp_end[index_i] = GET_DP_END(graph, w, index_i);
            beg_sn = dp_beg_sn[index_i] = (dp_beg[index_i]+pn)/pn-1; end_sn = dp_end_sn[index_i] = (dp_end[index_i]+pn)/pn-1;

            // loop query
            // init h, e                       
            _end_sn = MIN_OF_TWO(end_sn+1, dp_sn-1);
            for (i = beg_sn; i <= _end_sn; ++i) { // SIMD parallelization
                SIMDStorei(&dp_h[i], min_inf); SIMDStorei(&dp_e[i], min_inf); SIMDStorei(&dp_f[i], min_inf);
            }

            // get max m and e
            for (i = 0; i < pre_n[index_i]; ++i) {
                pre_i = pre_index[index_i][i];
                pre_dp_h = DP_HE + pre_i * 2 * dp_sn, pre_dp_e = pre_dp_h + dp_sn;
                pre_end = dp_end[pre_i];
                pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];

                // set M from (pre_i, q_i-1)
                _end_sn = MIN_OF_TWO(((pre_end + pn)/pn-1), dp_sn-1);
                SIMDi first = SIMDShiftRight(min_inf, SIMDTotalBytes-SIMDShiftOneNi16);
                for (sn_i = pre_beg_sn; sn_i <= _end_sn; ++sn_i) { // SIMD parallelization
                    a = SIMDLoadi(&pre_dp_h[sn_i]); b = SIMDLoadi(&dp_h[sn_i]);
                    SIMDi remain = SIMDShiftLeft(a, SIMDShiftOneNi16);
                    SIMDStorei(&dp_h[sn_i], SIMDMaxi16(SIMDOri(first, remain), b));
                    first = SIMDShiftRight(a, SIMDTotalBytes-SIMDShiftOneNi16);
                }

                // set E from (pre_i, q_i)
                for (sn_i = pre_beg_sn; sn_i <= pre_end_sn; ++sn_i) { // SIMD parallelization
                    a = SIMDLoadi(&pre_dp_e[sn_i]); b = SIMDLoadi(&dp_e[sn_i]);
                    SIMDStorei(&dp_e[sn_i], SIMDMaxi16(a, b));
                }
            }

            if (n_cigar && graph_cigar) {
                z = backtrack_z + (index_i-1) * d_sn;

                // compare M, E, and F
                for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { // SIMD parallelization
                    a = SIMDLoadi(&dp_h[sn_i]); b = SIMDLoadi(&q[sn_i]); a = SIMDAddi16(a, b);
                    c = SIMDLoadi(&dp_e[sn_i]); d = SIMDLoadi(&hd[sn_i]); 
                    SIMDGetIfGreateri16(d, a, c, a, he, hm);
                    SIMDStorei(&hd[sn_i], d); SIMDStorei(&dp_h[sn_i], a);
                }

                // set F from (index_i, q_i-1)
                //_dp_f[beg] = inf_min;
                for (q_i = beg+1; q_i <= end; ++q_i) { // no SIMD parallelization
                    _dp_f[q_i] = MAX_OF_TWO(_dp_h[q_i-1] - gap_OE, _dp_f[q_i-1] - gap_ext);
                }

                for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { // SIMD parallelization
                    a = SIMDLoadi(&dp_h[sn_i]); b = SIMDLoadi(&dp_f[sn_i]); c = SIMDLoadi(&dp_e[sn_i]); d = SIMDLoadi(&hd[sn_i]);
                    // h, hd for current cell
                    SIMDGetIfGreateri16(d, a, b, a, hf, d);
                    SIMDStorei(&dp_h[sn_i], a);

                    // e, ed for next cell
                    tmp_e = SIMDSubi16(a, gap_oe);
                    c = SIMDSubi16(c, gap_e);
                    if (gap_open) {
                        SIMDGetIfGreateri16(ed, c, c, tmp_e, ee, em);
                        SIMDStorei(&dp_e[sn_i], c);
                    } else {
                        SIMDSetIfEquali16(ed, d, hm, em, ee);
                        SIMDStorei(&dp_e[sn_i], tmp_e);
                    }

                    // fd for next cell
                    b = SIMDSubi16(b, gap_e);
                    if (gap_open) {
                        SIMDSetIfGreateri16(fd, b, tmp_e, ff, fm);
                    } else {
                        SIMDSetIfEquali16(fd, d, hm, fm ,ff);
                    }
                    SIMDStorei(&z[sn_i], SIMDOri(SIMDOri(d, ed), fd));
                }
            } else {
                // compare M, E, and F
                for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { // SIMD parallelization
                    a = SIMDLoadi(&dp_h[sn_i]); b = SIMDLoadi(&q[sn_i]); c = SIMDLoadi(&dp_e[sn_i]);
                    a = SIMDAddi16(a, b);  // q[0] == 0
                    SIMDStorei(&dp_h[sn_i], SIMDMaxi16(a, c));
                    print_simd(&dp_h[sn_i], "dp_h", int16_t);
                }

                // set F from (index_i, q_i-1)
                _dp_f[beg] = inf_min;
                for (q_i = beg+1; q_i <= end; ++q_i) { // no SIMD parallelization
                    _dp_f[q_i] = MAX_OF_TWO(_dp_h[q_i-1] - gap_OE, _dp_f[q_i-1] - gap_ext);
                }

                for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { // SIMD parallelization
                    a = SIMDLoadi(&dp_h[sn_i]); b = SIMDLoadi(&dp_f[sn_i]); c = SIMDLoadi(&dp_e[sn_i]);
                    // h for current cell
                    a = SIMDMaxi16(a,b);
                    SIMDStorei(&dp_h[sn_i], a);

                    // e for next cell
                    tmp_e = SIMDSubi16(a, gap_oe);
                    c = SIMDSubi16(c, gap_e);
                    SIMDStorei(&dp_e[sn_i], SIMDMaxi16(c, tmp_e));
                }
            }
        }
    }

#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("index: %d\t", j);
        dp_h = DP_HE + j * 2 * dp_sn, dp_e = dp_h + dp_sn;
        _dp_h = (int16_t*)dp_h, _dp_e = (int16_t*)dp_e;
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {
        //for (i = 0; i <= qlen; ++i) {
            printf("%d:(%d,%d)\t", i, _dp_h[i], _dp_e[i]);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    for (i = 0; i < graph->node[POA_SINK_NODE_ID].in_edge_n; ++i) { // for global alignment, find best backtrack position
        in_id = graph->node[POA_SINK_NODE_ID].in_id[i];
        in_index = bpoa_graph_node_id_to_index(graph, in_id);
        dp_h = DP_HE + in_index * 2 * dp_sn;
        _dp_h = (int16_t*)dp_h;
        _set_max_score(best_score, best_i, best_j, _dp_h[qlen], in_index, qlen);
    }

    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
    { // backtrack from best score
        if (n_cigar && graph_cigar) {
            simd_bpoa_backtrack(int16_t, DP_HE, dp_sn, (int16_t)bpt->match, (int16_t)bpt->mismatch, (int16_t)bpt->gap_ext, pre_index, pre_n, backtrack_z, d_sn, best_i, best_j, graph, query, n_cigar, graph_cigar);
        }
    }

    { // free variables
        SIMDFree(qp); free(dp_beg); free(dp_end); free(dp_beg_sn); free(dp_end_sn);
        if (n_cigar && graph_cigar) SIMDFree(hd);
        for (i = 0; i < graph->node_n; ++i) free(pre_index[i]); 
        free(pre_index); free(pre_n);
    }
    return best_score;
}
 
// based on min,max and w
// set score_n and id_n (8, 16, 32 bits)
int simd_bpoa_banded_global_align_sequence_with_graph(bpoa_graph_t *graph, uint8_t *query, int qlen, bpoa_para_t *bpt, int *n_cigar, bpoa_cigar_t **graph_cigar) {
    //bpt->simd_flag = 0;
    if (bpt->simd_flag == 0) return bpoa_banded_global_align_sequence_with_graph(graph, query, qlen, bpt, n_cigar, graph_cigar);

    int max_score;
    if (bpt->simd_flag & SIMD_AVX512F && !(bpt->simd_flag & SIMD_AVX512BW)) max_score = 32768; // AVX512F has no 8/16 bits operations
    else {
        int max_l = bpoa_graph_node_id_to_max_remain(graph, POA_SRC_NODE_ID);
        int min_l = bpoa_graph_node_id_to_min_remain(graph, POA_SRC_NODE_ID);
        if (min_l <= qlen && qlen <= max_l) {
            max_score = qlen * bpt->m;
        } else if (max_l < qlen) {
            max_score = qlen * bpt->m - bpt->gap_open - (max_l-qlen) * bpt->gap_ext;
        } else max_score = qlen * bpt->m - bpt->gap_open - (qlen-min_l) * bpt->gap_ext;
        max_score = MAX_OF_TWO(max_score, bpt->gap_open + max_l * bpt->gap_ext + bpt->gap_open + qlen * bpt->gap_ext);
    }
    int _best_score;
    if (max_score <= 127 - bpt->mismatch - bpt->gap_open - bpt->gap_ext) { // DP_H/E/F: 8  bits
        _simd_p8.inf_min = MAX_OF_TWO(INF_8_MIN + bpt->mismatch, INF_8_MIN + bpt->gap_open + bpt->gap_ext);
        simd_bpoa_banded_global_align_sequence_with_graph_core(graph, query, qlen, bpt, n_cigar, graph_cigar, _simd_p8, int8_t, SIMDSetOnei8, SIMDMaxi8, SIMDAddi8, SIMDSubi8, SIMDShiftOneNi8, SIMDSetIfGreateri8, SIMDGetIfGreateri8, SIMDSetIfEquali8, _best_score);
        return _best_score;
        //return simd8_bpoa_banded_global_align_sequence_with_graph_core(graph, query, qlen, bpt, n_cigar, graph_cigar, _simd_p8);
    } else if (max_score <= 32767 - bpt->mismatch - bpt->gap_open - bpt->gap_ext) { // DP_H/E/F: 16 bits
        _simd_p16.inf_min = MAX_OF_TWO(INF_16_MIN + bpt->mismatch, INF_16_MIN + bpt->gap_open + bpt->gap_ext);
        simd_bpoa_banded_global_align_sequence_with_graph_core(graph, query, qlen, bpt, n_cigar, graph_cigar, _simd_p16, int16_t, SIMDSetOnei16, SIMDMaxi16, SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, SIMDSetIfEquali16, _best_score);
        return _best_score;
        //return simd16_bpoa_banded_global_align_sequence_with_graph_core(graph, query, qlen, bpt, n_cigar, graph_cigar, _simd_p16);
    } else { // 2147483647, DP_H/E/F: 32 bits
        _simd_p32.inf_min = MAX_OF_TWO(INF_32_MIN + bpt->mismatch, INF_32_MIN + bpt->gap_open + bpt->gap_ext);
        simd_bpoa_banded_global_align_sequence_with_graph_core(graph, query, qlen, bpt, n_cigar, graph_cigar, _simd_p32, int32_t, SIMDSetOnei32, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, SIMDSetIfEquali32, _best_score);
        return _best_score;
        //return simd32_bpoa_banded_global_align_sequence_with_graph_core(graph, query, qlen, bpt, n_cigar, graph_cigar, _simd_p32);
    }
    return 0;
}
