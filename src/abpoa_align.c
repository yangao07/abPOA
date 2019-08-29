#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "abpoa_align.h"
#include "utils.h"

#define TODO_M   0x1
#define TODO_X   0x2
#define TODO_MX  0x3
#define TODO_E   0x4
#define TODO_ME  0x5
#define TODO_F   0x8
#define TODO_MF  0x9
#define TODO_MEF 0xd

#define HAS_M    0x10
#define HAS_X    0x20
#define HAS_MX   0x30
#define HAS_E    0x40
#define HAS_ME   0x50
#define HAS_F    0x80


// xl, el, fl: current extension length of X, E and F
// h: max, e: vertical, f: horizontal
typedef struct {
    int64_t h:32, xl:32;
    int64_t e:32, el:32;
    int64_t f:32, fl:32; 
    uint8_t todo_map; // XXX need to stored seperately
} dp_matrix_t;

void gen_simple_mat(int m, int *mat, int match, int mismatch) {
    int i, j;
    match = match < 0 ? -match : match;
    mismatch = mismatch > 0? -mismatch : mismatch;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j ? match : mismatch;
        mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = 0;
}

void abpoa_set_gap_mode(abpoa_para_t *abpt) {
    if (abpt->gap_open1 == 0) abpt->gap_mode = ABPOA_LINEAR_GAP;
    else if (abpt->gap_open1 > 0 && abpt->gap_open2 == 0) abpt->gap_mode = ABPOA_AFFINE_GAP;
    else abpt->gap_mode = ABPOA_CONVEX_GAP;
}

abpoa_para_t *abpoa_init_para(void) {
    abpoa_para_t *abpt = (abpoa_para_t*)_err_malloc(sizeof(abpoa_para_t));
    abpt->align_mode = ABPOA_GLOBAL_MODE;
    abpt->gap_mode = ABPOA_CONVEX_GAP;
    abpt->zdrop = -1;     // disable zdrop
    abpt->end_bonus = -1; // disable end bouns
    abpt->bw = -1;        // disable bandwidth
    abpt->use_ada = 0;    // use adaptive band
    abpt->ret_cigar = 1;  // return cigar
    abpt->out_msa = 0;    // output msa
    abpt->out_cons = 0;   // output consensus sequence in msa
    abpt->multip = ABPOA_MULTIP; // muliple consensus for multiploid data
    abpt->min_fre = 0.3;
    abpt->cons_agrm = ABPOA_HB;   // consensus calling algorithm 
    abpt->use_read_ids = 0;
    abpt->out_pog= 0; // generate partial order graph

    // number of residue types
    abpt->m = 5; // nucleotides
    abpt->mat = NULL; // TODO score matrix for aa

    // score matrix
    abpt->match = ABPOA_MATCH;
    abpt->mismatch = ABPOA_MISMATCH;
    abpt->gap_open1 = ABPOA_GAP_OPEN1;
    abpt->gap_open2 = ABPOA_GAP_OPEN2;
    abpt->gap_ext1 = ABPOA_GAP_EXT1;
    abpt->gap_ext2 = ABPOA_GAP_EXT2;

    abpt->simd_flag = simd_check();

    return abpt;
}

void abpoa_free_para(abpoa_para_t *abpt) {
    if (abpt->mat != NULL) free(abpt->mat);
    free(abpt);
}

// 0. init score table
// 1. traverse graph with topological order; banded width
// 2. backtrack and generate graph_cigar
// h: max score, e: score from vertical, f: score from horizontal
int abpoa_global_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar) {
    dp_matrix_t *dp_matrix; // Full: (tlen + 1) * (qlen + 1); Banded: (tlen+1) * (2 * w + 1)
    int *qp, *mat = abpt->mat; // query profile
    int i, j, k, gap_o = abpt->gap_open1, gap_e = abpt->gap_ext1, gap_oe = abpt->gap_open1 + abpt->gap_ext1, w, matrix_row_n, matrix_col_n, z_col_n;
    int node_id, target_node_n = graph->node_n - 2; // exclude start and end nodes
    uint64_t *z; // backtrack matrix; in each cell: hd << 33 | ed << 2 | fd
    //                                 h<<62|h_id<<33|e<<31|e_i<<2|f
    //                                 h_id/e_i: 29 bit XXX cause error when in_edge_n >= pow(2,29)
    //                                 MATCH:0, DELETION:1, INSERTION:2, MISMATCH:3
    int best_score = INT32_MIN, best_i=0, best_j=0;

    // allocate memory 
    matrix_row_n = target_node_n + 1, matrix_col_n = qlen + 1;
    dp_matrix = (dp_matrix_t*)_err_calloc(matrix_row_n * matrix_col_n, sizeof(dp_matrix_t));
    w = abpt->bw <= 0 ? qlen : abpt->bw; // band width TODO adaptive band width, shift with diagonal
    z_col_n = qlen < 2 * w + 1 ? qlen : 2 * w + 1; // TODO adaptive band width
    z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0;
    qp = (int*)_err_malloc(qlen * abpt->m);

    // generate the query profile
    for (k = i = 0; k < abpt->m; ++k) {
        const int *p = &mat[k * abpt->m];
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }
    // fill the first row
    dp_matrix[0].h = 0; dp_matrix[0].e = INT32_MIN; // f is useless
    for (j = 1; j < matrix_col_n && j <= w; ++j) 
        dp_matrix[j].h = dp_matrix[j].f = -(gap_o + gap_e * j), dp_matrix[j].e = INT32_MIN;
    for (; j < matrix_col_n; ++j) dp_matrix[j].h = dp_matrix[j].e = INT32_MIN; // everything outside the band is -INF
    // DP loop
    for (i = 1; i <= target_node_n; ++i) { // target graph is in the outer loop
        node_id = abpoa_graph_index_to_node_id(graph, i); // i: node index, node_i: node id
        int32_t f = INT32_MIN, beg, beg_, end;
        int *q = &qp[graph->node[node_id].base * qlen];
        beg = i > w ? i - w : 0; // set band boundary FIXME for multi-in node
        end = i + w + 1 < qlen + 1 ? i + w + 1: qlen + 1; // only loop through [beg,end) of the query sequence 
        int pre_n = graph->node[node_id].in_edge_n;
        int *pre_i = graph->node[node_id].in_id;
        int *pre_index = (int*)_err_malloc(pre_n * sizeof(int));
        for (j = 0; j < pre_n; ++j)
            pre_index[j] = abpoa_graph_node_id_to_index(graph, pre_i[j]);

        dp_matrix_t *pre_matrix;
        dp_matrix_t *cur_matrix = &dp_matrix[i * matrix_col_n];

        if (n_cigar && graph_cigar) { // store backtrack direction
            uint64_t *zi = &z[(long)(i-1) * z_col_n];
            int32_t m, h, e, tmp; uint64_t pre_i, m_i, fd, ed, e_i, hd, mx, m0=0, e1=1, f2=2, x3=3; // h<<62|h_id<<33|e<<31|e_i<<2|f
            // fill the first column
            if (beg == 0) {
                int min_rank = abpoa_graph_node_id_to_min_rank(graph, node_id);
                cur_matrix[beg].h = -(gap_o + gap_e * min_rank);
                cur_matrix[beg].e = -(gap_o + gap_e * min_rank);
                cur_matrix[beg].f = INT32_MIN;
                beg_ = 1;
            } else beg_ = beg;
            for (j = beg_; j < end; ++j) {
                k = 0; // first precursor
                pre_matrix = &dp_matrix[pre_index[k] * matrix_col_n];
                m_i = e_i = pre_index[k];
                //cur_matrix[j].f = f = MAX_OF_TWO(cur_matrix[j-1].f, cur_matrix[j-1].h-gap_o) - gap_e;
                tmp = cur_matrix[j-1].h - gap_oe;
                f = cur_matrix[j-1].f - gap_e;
                cur_matrix[j].f = f > tmp ? f : tmp;
                f = cur_matrix[j].f;

                //cur_matrix[j].e = e = MAX_OF_TWO(pre_matrix[j].e, pre_matrix[j].h-gap_o) - gap_e;
                tmp = pre_matrix[j].h - gap_oe;
                e = pre_matrix[j].e - gap_e;
                cur_matrix[j].e = e > tmp ? e : tmp;
                e = cur_matrix[j].e;

                // XXX does not select best_M, best_E first. so does not prefer to M than E or F
                m = pre_matrix[j-1].h + q[j-1];
                mx = q[j-1] == abpt->match ? m0 : x3;
                // hd for current cell
                hd = m >= e ? (mx << 62 | m_i << 33) : (e1 << 62 | e_i << 33);
                h = m >= e ? m : e;
                hd = h >= f ? hd : f2 << 62 ;
                h = h >= f ? h : f;
                cur_matrix[j].h = h;

                // ed and fd for next cell
                ed = e - gap_e > m - gap_oe ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
                fd = f - gap_e > m - gap_oe ? f2 : 0;
                // if has multiple precursors
                for (k = 1; k < pre_n; ++k) {
                    pre_matrix = &dp_matrix[pre_index[k] * matrix_col_n];
                    pre_i = pre_index[k];

                    e = pre_matrix[j].h - gap_oe;
                    if (e > cur_matrix[j].e) {
                        cur_matrix[j].e = e;
                        e_i = pre_i;
                    }
                    e = pre_matrix[j].e - gap_e;
                    if (e > cur_matrix[j].e) {
                        cur_matrix[j].e = e;
                        e_i = pre_i;
                    }
                    e = cur_matrix[j].e;

                    if (pre_matrix[j-1].h + q[j-1] > m) {
                        m = pre_matrix[j-1].h + q[j-1];
                        m_i = pre_i;
                        mx = q[j-1] == abpt->match ? m0 : x3;
                    }
                    // hd for current cell
                    hd = m > cur_matrix[j].h ? (mx << 62 | m_i << 33) : hd;
                    cur_matrix[j].h = m > cur_matrix[j].h ? m : cur_matrix[j].h;

                    hd = cur_matrix[j].h >= e ? hd : (e1 << 62 | e_i << 33);
                    cur_matrix[j].h = cur_matrix[j].h >= e ? cur_matrix[j].h : e;

                    // ed and fd for next cell
                    ed = e - gap_e > m - gap_oe ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
                    fd = f - gap_e > m - gap_oe ? f2 : 0;
                }
                zi[j - beg_] = hd | ed | fd;
                // printf("direction(%d, %d, %d): hd: %d, h_id: %d, ed: %d, e_id: %d, f_id: %d\n", i, j, beg_, (int)(hd >> 62), (int)((hd >> 33) & 0x1fffffff), (int)(ed >> 31), (int)((ed >> 2) & 0x1fffffff), (int)fd);
            }
            // fill the last column for next row
            if (end <= qlen) {
                cur_matrix[end].e = INT32_MIN;
                cur_matrix[end].h = cur_matrix[end].f = MAX_OF_TWO(cur_matrix[end-1].e, cur_matrix[end-1].h) - gap_e;
            }
        } else { // only compute score
            int32_t m, h, e;
            // fill the first column
            if (beg == 0) {
                cur_matrix[beg].h = -(gap_o + gap_e * i);
                cur_matrix[beg].e = -(gap_o + gap_e * i);
                cur_matrix[beg].f = INT32_MIN;
                beg_ = 1;
            } else beg_ = beg;
            for (j = beg_; j < end; ++j) {
                k = 0; // first precursor
                pre_matrix = &dp_matrix[pre_index[k] * matrix_col_n];
                m = pre_matrix[j-1].h + q[j-1];
                cur_matrix[j].e = e = MAX_OF_TWO(pre_matrix[j].e, pre_matrix[j].h-gap_o) - gap_e;
                cur_matrix[j].f = f = MAX_OF_TWO(cur_matrix[j-1].f, cur_matrix[j-1].h-gap_o) - gap_e;
                h = MAX_OF_THREE(m, e, f);
                cur_matrix[j].h = h;
                for (k = 1; k < pre_n; ++k) {
                    pre_matrix = &dp_matrix[pre_index[k] * matrix_col_n];
                    m = pre_matrix[j-1].h + q[j-1];
                    e = MAX_OF_TWO(cur_matrix[j].e, MAX_OF_TWO(pre_matrix[j].e, pre_matrix[j].h-gap_o) - gap_e);
                    f = MAX_OF_TWO(cur_matrix[j].f, MAX_OF_TWO(cur_matrix[j-1].f, cur_matrix[j-1].h-gap_o) - gap_e);
                    h = MAX_OF_THREE(m, e, f);

                    cur_matrix[j].e = e;
                    cur_matrix[j].f = f;
                    cur_matrix[j].h = h;
                }
            }
            // fill the last column for next row
            if (end <= qlen) {
                cur_matrix[end].e = INT32_MIN;
                cur_matrix[end].h = cur_matrix[end].f = MAX_OF_TWO(cur_matrix[end-1].e, cur_matrix[end-1].h) - gap_e;
            }
        }
        free(pre_index);
    }
#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("index: %d\t", j);
        for (i = 0; i <= qlen; ++i) {
            printf("%d:(%d,%d,%d)\t", i, dp_matrix[j*matrix_col_n+i].h, dp_matrix[j*matrix_col_n+i].e, dp_matrix[j*matrix_col_n+i].f);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    for (i = 0; i < graph->node[ABPOA_SINK_NODE_ID].in_edge_n; ++i) {
        in_id = graph->node[ABPOA_SINK_NODE_ID].in_id[i];
        in_index = abpoa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, dp_matrix[in_index * matrix_col_n + qlen].h, in_index, qlen);
    }

#ifdef __DEBUG__
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
#endif
    // backtrack from best score
    if (n_cigar && graph_cigar) {
        int n_c = 0, m_c = 0, id, which;
        int op_shift[4] = {62, 31, 0, 62}, id_shift[4] = {33, 2, 0, 33};
        uint64_t d;
        abpoa_cigar_t *cigar = 0;
        i = best_i, j = best_j, id = abpoa_graph_index_to_node_id(graph, i), which = 0;
        while (i > 0 && j > 0) {
            d = z[(long)(i-1) * z_col_n + (j-1 - (i-1 > w ? i-1 - w : 0))];
            which = (d >> op_shift[which]) & 3;
            // printf("(%d,%d) %c\n", i, j, ABPOA_CIGAR_STR[which]);
            if (which == 0) { // match
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id/*index_to_id*/, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = abpoa_graph_index_to_node_id(graph, i);
                j--;
            } else if (which == 3) { // mismatch
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDIFF, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = abpoa_graph_index_to_node_id(graph, i);
                j--;
            } else if (which == 1) { // deletion
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = abpoa_graph_index_to_node_id(graph, i);
            } else { // insertion
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1);
                j--;
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
    free(dp_matrix); free(qp); free(z);
    return best_score;
}

// TODO extend
int abpoa_ada_extend_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar) {
    dp_matrix_t *dp_matrix, *cur_line, *pre_line, *next_line, *cur_dp; // Full: (tlen + 1) * (qlen + 1); Banded: (tlen+1) * (2 * w + 1)
    int *qp, *mat = abpt->mat, gap_o = abpt->gap_open1, gap_e = abpt->gap_ext1, gap_oe = abpt->gap_open1 + abpt->gap_ext1; // query profile
    int **pre_index, *pre_n, pre_i, **next_index, *next_n, next_i;
    int node_id, index_i, min_rank, q_i;
    int i, j, k, w, *dp_beg_cen, *dp_end_cen;
    int target_node_n = graph->node_n - 2; // exclude start and end nodes
    int matrix_row_n = graph->node_n, matrix_col_n = qlen + 2, z_col_n = qlen; // TODO use less memory ???
    uint64_t *backtrack_z, *z; // backtrack matrix; in each cell: hd << 33 | ed << 2 | fd
    //                                 h<<62|h_id<<33|e<<31|e_i<<2|f
    //                                 h_id/e_i: 29 bit XXX cause error when in_edge_n >= pow(2,29)
    //                                 MATCH:0, DELETION:1, INSERTION:2, MISMATCH:3
    int best_score = INT32_MIN, best_i=0, best_j=0;

    // allocate memory 
    dp_matrix = (dp_matrix_t*)_err_calloc(matrix_row_n * matrix_col_n, sizeof(dp_matrix_t)); // make sure no invalid write 
    dp_beg_cen = (int*)_err_malloc((matrix_row_n) * sizeof(int)); for (i = 1; i < matrix_row_n; ++i) dp_beg_cen[i] = qlen;
    dp_end_cen = (int*)_err_calloc(matrix_row_n, sizeof(int));
    // TODO if no band is used
    // w = abpt->bw <= 0 ? (MIN_OF_TWO(graph->rank_n-2, qlen)) : (MIN_OF_THREE(graph->rank_n-2, qlen, abpt->bw));
    w = 300; // XXX small may cause NO alignment result XXX
    backtrack_z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0;

    qp = (int*)_err_malloc(qlen * abpt->m); // TODO score profile for diagnal

    pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));
    next_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    next_n = (int*)_err_malloc(graph->node_n * sizeof(int*));

    // generate the query profile
    for (k = i = 0; k < abpt->m; ++k) {
        const int *p = &mat[k * abpt->m];
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }

    // index of pre-node
    for (i = 0; i < graph->node_n; ++i) {
        node_id = abpoa_graph_index_to_node_id(graph, i); // i: node index
        pre_n[i] = graph->node[node_id].in_edge_n;
        pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
        for (j = 0; j < pre_n[i]; ++j) {
            pre_index[i][j] = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
        }
        
        next_n[i] = graph->node[node_id].out_edge_n;
        next_index[i] = (int*)_err_malloc(next_n[i] * sizeof(int));
        for (j = 0; j < next_n[i]; ++j) {
            next_index[i][j] = abpoa_graph_node_id_to_index(graph, graph->node[node_id].out_id[j]);
        }
    }

    // TODO no need to store F ???
    // cur_cell: H[i,j], E[i+1,j], F[i,j+1]
    // fill the first row
    cur_line = dp_matrix;
    cur_line->h = 0; cur_line->e = -gap_oe; cur_line->f = -gap_oe;
    cur_line->xl = 0; cur_line->fl = cur_line->el = 1; 
    cur_line->todo_map = HAS_MX; 
    for (i = 1; i <= w; ++i) {
        cur_dp = cur_line + i;
        cur_dp->e = INT32_MIN; cur_dp->h = -gap_o - gap_e * i; cur_dp->f = cur_dp->h - gap_e;
        cur_dp->fl = i+1; cur_dp->el = cur_dp->xl = 0; 
        cur_dp->todo_map = HAS_MX; 
    }
    for (i = 0; i < next_n[ABPOA_SRC_NODE_ID]; ++i) {
        next_i = next_index[ABPOA_SRC_NODE_ID][i];
        next_line = dp_matrix + next_i * matrix_col_n;
        next_line[0].todo_map = TODO_E;
        for (j = 1; j<= w+1; ++j) {
            next_line[j].todo_map = TODO_MF;
        }
        dp_beg_cen[next_i] = 1;
        dp_end_cen[next_i] = 1;
    }

    // DP loop
    int32_t tmp, h, m, e, f, beg, end, _beg, _end, _next_cen;
    uint64_t m_i, fd, ed, e_i, hd, mx, m0=0x0, e1=0x1, f2=0x2, x3=0x3; // h<<62|h_id<<33|e<<31|e_i<<2|f TODO f<<62|e<<60|e_id<<31|h<<29|h_id
                                                                       // m0=0x0 << 62, x3=0x3<<62
    dp_matrix_t *m_pre_dp, *e_pre_dp, *f_pre_dp;

    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {
        cur_line = &dp_matrix[index_i * matrix_col_n];
        node_id = abpoa_graph_index_to_node_id(graph, index_i);
        min_rank = abpoa_graph_node_id_to_min_rank(graph, node_id);
        int *q = &qp[graph->node[node_id].base * qlen];
        beg = dp_beg_cen[index_i] >= w ? dp_beg_cen[index_i] - w : 0;
        end = dp_end_cen[index_i] + w <= qlen ? dp_end_cen[index_i] + w : qlen;

        if (beg == 0) {
            cur_line[0].h = -(gap_o + gap_e * min_rank);
            cur_line[0].e = -(gap_o + gap_e * (min_rank+1));
            cur_line[0].f = INT32_MIN;
            _beg = 1;
        } else _beg = beg;
        for (q_i = _beg; q_i <= end; ++q_i) {
            if (!cur_line[q_i].todo_map) continue;
            cur_dp = cur_line + q_i;
            m = e = f = INT32_MIN;
            m_i = e_i = -1;
            z = &backtrack_z[(index_i-1) * qlen];

            if (cur_dp->todo_map & TODO_M) { // M: from (pre_i, q_i-1)
                for (i = 0; i < pre_n[index_i]; ++i) {
                    pre_i = pre_index[index_i][i];
                    pre_line = dp_matrix + pre_i * matrix_col_n;
                    m_pre_dp = pre_line + q_i-1;
                    if ((m_pre_dp->todo_map & HAS_M) && m_pre_dp->h > m) {
                        m = m_pre_dp->h; m_i = pre_i;
                    }
                }
            }
            if (cur_dp->todo_map & TODO_E) { // E: from (pre_i, q_i)
                for (j = 0; j < pre_n[index_i]; ++j) {
                    pre_i = pre_index[index_i][j];
                    pre_line = dp_matrix + pre_i * matrix_col_n;
                    e_pre_dp = pre_line + q_i;
                    if ((e_pre_dp->todo_map & HAS_E) && e_pre_dp->e > e) {
                        e = e_pre_dp->e; e_i = pre_i;
                    }
                }
            }
            if (cur_dp->todo_map & TODO_F) { // F: from (index_i, q_i-1)
                f_pre_dp = cur_line + q_i - 1;
                f = f_pre_dp->f;
            }

            m += q[q_i-1]; mx = q[q_i-1] == abpt->match ? m0 : x3;
            // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell
            // TODO z-drop
             
            // h,hd for current cell
            hd = m >= e ? (mx << 62 | m_i << 33) : (e1 << 62 | e_i << 33);
            h = m >= e ? m : e;
            hd = h >=f ? hd : f2 << 62;
            h = h >= f ? h : f;
            cur_dp->h = h;

            // e,ed for next cell
            tmp = h - gap_oe;
            e -= gap_e;
            ed = e > tmp ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
            e = e > tmp ? e : tmp;
            cur_dp->e = e;

            // f, fd for next cell
            f -= gap_e;
            fd = f > tmp ? f2 : 0;
            f = f > tmp ? f : tmp;
            cur_dp->f = f;

            z[q_i-1] = hd | ed | fd;
        }
        // set TODO_M/E/F for next_line based one cur_line[start] and cur_line[end]
        // set start, end
        if (cur_line[beg].h >= cur_line[end].h) { // move downward
            _next_cen = dp_beg_cen[index_i];
            _beg = beg;
            _end = dp_beg_cen[index_i] + w <= qlen ? dp_beg_cen[index_i] + w : qlen;
            for (i = 0; i < next_n[index_i]; ++i) {
                next_i = next_index[index_i][i];
                next_line = dp_matrix + next_i * matrix_col_n;
                dp_beg_cen[next_i] = dp_beg_cen[next_i] < _next_cen ? dp_beg_cen[next_i] : _next_cen;
                dp_end_cen[next_i] = dp_end_cen[next_i] > _next_cen ? dp_end_cen[next_i] : _next_cen;

                for (j = beg; j <= _end-1; ++j)
                    cur_line[j].todo_map |= HAS_ME;
                cur_line[_end].todo_map |= HAS_E;

                next_line[_beg].todo_map |= TODO_E;
                for (j = _beg+1; j <= _end; ++j) {
                    next_line[j].todo_map |= TODO_MEF;
                }
            }
        } else { // move rightward
            _next_cen = dp_end_cen[index_i] >= qlen ? qlen : dp_end_cen[index_i] + 1;
            _beg = _next_cen >= w ? _next_cen - w : 0;
            _end = end < qlen ? end + 1 : qlen;
            for (i = 0; i < next_n[index_i]; ++i) {
                next_i = next_index[index_i][i];
                next_line = dp_matrix + next_i * matrix_col_n;
                dp_end_cen[next_i] = dp_end_cen[next_i] > _next_cen ? dp_end_cen[next_i] : _next_cen;
                dp_beg_cen[next_i] = dp_beg_cen[next_i] < _next_cen ? dp_beg_cen[next_i] : _next_cen;

                cur_line[_beg-1].todo_map |= HAS_M;
                for (j = _beg; j <= _end-1; ++j) {
                    cur_line[j].todo_map |= HAS_ME;
                }
                if (end == qlen) cur_line[end].todo_map |= HAS_E;

                next_line[_beg].todo_map |= TODO_ME;
                for (j = _beg+1; j < _end; ++j) {
                    next_line[j].todo_map |= TODO_MEF;
                }
                next_line[_end].todo_map |= end == qlen ? TODO_MEF : TODO_MF;
            }
        }
    }

#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("rank: %d\t", j);
        for (i = 0; i <= qlen; ++i) {
            printf("%d:(%d,%d,%d)\t", i, dp_matrix[j*matrix_col_n+i].h, dp_matrix[j*matrix_col_n+i].e, dp_matrix[j*matrix_col_n+i].f);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    // for global alignment, find best backtrack position
    // TODO semi-global, local
    for (i = 0; i < graph->node[ABPOA_SINK_NODE_ID].in_edge_n; ++i) {
        in_id = graph->node[ABPOA_SINK_NODE_ID].in_id[i];
        in_index = abpoa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, dp_matrix[in_index * matrix_col_n + qlen].h, in_index, qlen);
    }

#ifdef __DEBUG__
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
#endif
    // backtrack from best score
    if (n_cigar && graph_cigar) {
        int n_c = 0, m_c = 0, id, which;
        int op_shift[4] = {62, 31, 0, 62}, id_shift[4] = {33, 2, 0, 33};
        uint64_t d;
        abpoa_cigar_t *cigar = 0;
        i = best_i, j = best_j, id = abpoa_graph_index_to_node_id(graph, i), which = 0;
        while (i > 0 && j > 0) {
            d = backtrack_z[(long)(i-1) * z_col_n + j-1];
            which = (d >> op_shift[which]) & 3;
            // printf("(%d,%d) %c\n", i, j, ABPOA_CIGAR_STR[which]);
            if (which == 0) { // match
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = abpoa_graph_index_to_node_id(graph, i);
                j--;
            } else if (which == 3) { // mismatch
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDIFF, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = abpoa_graph_index_to_node_id(graph, i);
                j--;
            } else if (which == 1) { // deletion
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = abpoa_graph_index_to_node_id(graph, i);
            } else { // insertion
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1);
                j--;
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
    free(dp_matrix); free(qp); free(backtrack_z); free(dp_beg_cen); free(dp_end_cen);
    for (i = 0; i < graph->node_n; ++i) { free(pre_index[i]); } free(pre_index); free(pre_n);
    for (i = 0; i < graph->node_n; ++i) { free(next_index[i]); } free(next_index); free(next_n);
    return best_score;

    return 0;
}

/*int abpoa_banded_global_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar) {
    int8_t *qp, *mat = abpt->mat, gap_o = abpt->gap_open1, gap_e = abpt->gap_ext1, gap_oe = abpt->gap_open1 + abpt->gap_ext1; // query profile
    int **pre_index, *pre_n, pre_i, **next_index, *next_n;
    int node_id, index_i, q_i;
    int i, j, k, w, *dp_beg, *dp_end;
    int target_node_n = graph->node_n - 2; // exclude start and end nodes
    int matrix_row_n = graph->node_n, matrix_col_n = qlen + 1, z_col_n = qlen; // TODO use less memory? doable, casue band range is known
    int *DP_H, *DP_E, *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *dp_f;

    uint64_t *backtrack_z, *z; // backtrack matrix; in each cell: hd << 33 | ed << 2 | fd
    uint64_t *hd, *mx, *m_pre_i, *e_pre_i;
    //                                 h<<62|h_id<<33|e<<31|e_i<<2|f
    //                                 h_id/e_i: 29 bit XXX cause error when in_edge_n >= pow(2,29)
    //                                 MATCH:0, DELETION:1, INSERTION:2, MISMATCH:3
    int best_score = INT32_MIN, best_i=0, best_j=0;

    // allocate memory 
    DP_H = (int*)_err_malloc(matrix_row_n * matrix_col_n * sizeof(int));
    DP_E = (int*)_err_malloc(matrix_row_n * matrix_col_n * sizeof(int));
    dp_f = (int*)_err_malloc(matrix_col_n * sizeof(int));

    dp_beg = (int*)_err_malloc((matrix_row_n) * sizeof(int)); dp_end = (int*)_err_calloc(matrix_row_n, sizeof(int));
    // if w <= 0, do whole global
    // w = abpt->bw <= 0 ? qlen;
    w = 50;
    // calculate band range for each row:
    // have: min_rank, max_rank, min_remain, max_remain
    // then: min_len = min_rank + min_remain, max_len = min_rank + max_remain
    // range: (min_of_two(min_rank, min_rank+qlen-max_len), max_of_two(min_rank+qlen-min_len, max_rank))
    // with w: (min-w, max+w)

    backtrack_z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0;
    qp = (int8_t*)_err_malloc(qlen * abpt->m);

    hd = (uint64_t*)_err_malloc(qlen * sizeof(int64_t));
    mx = (uint64_t*)_err_malloc(qlen * sizeof(int64_t));
    m_pre_i = (uint64_t*)_err_malloc(matrix_col_n * sizeof(int64_t));
    e_pre_i = (uint64_t*)_err_malloc(matrix_col_n * sizeof(int64_t));

    pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));
    next_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    next_n = (int*)_err_malloc(graph->node_n * sizeof(int*));

    // generate the query profile
    for (k = i = 0; k < abpt->m; ++k) {
        const int8_t *p = &mat[k * abpt->m];
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }

    // index of pre-node
    for (i = 0; i < graph->node_n; ++i) {
        node_id = abpoa_graph_index_to_node_id(graph, i); // i: node index
        pre_n[i] = graph->node[node_id].in_edge_n;
        pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
        for (j = 0; j < pre_n[i]; ++j) {
            pre_index[i][j] = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
        }
        
        next_n[i] = graph->node[node_id].out_edge_n;
        next_index[i] = (int*)_err_malloc(next_n[i] * sizeof(int));
        for (j = 0; j < next_n[i]; ++j) {
            next_index[i][j] = abpoa_graph_node_id_to_index(graph, graph->node[node_id].out_id[j]);
        }
    }

    // DP cell: H[i,j], E[i+1,j], F[i,j+1]
    // fill the first row
    dp_beg[0] = GET_DP_BEGIN(graph, w, 0); dp_end[0] = GET_DP_END(graph, w, 0);

    dp_h = DP_H; dp_e = DP_E;
    dp_h[0] = 0; dp_e[0] = -gap_oe;
    for (i = 1; i <= dp_end[0]; ++i) {
        dp_e[i] = INT32_MIN; dp_h[i] = -gap_o - gap_e * i;
    }

    // DP loop
    int tmp, beg, end, _beg, _end, pre_beg, pre_end;
    uint64_t fd, ed; // f|e|e_id|h|h_id
    uint64_t m0=0x0, e1=0x1, f2=0x2, x3=0x3, he, hf, ee, ff;
                                                     
    he = e1 << HOP_OFF_SET, hf = f2 << HOP_OFF_SET;
    ee = e1 << EOP_OFF_SET; ff = f2 << FOP_OFF_SET;

    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {
        node_id = abpoa_graph_index_to_node_id(graph, index_i);

        int8_t *q = &qp[graph->node[node_id].base * qlen];
        dp_h = DP_H + index_i * matrix_col_n; dp_e = DP_E + index_i * matrix_col_n;
        z = &backtrack_z[(index_i-1) * qlen];

        dp_beg[index_i] = GET_DP_BEGIN(graph, w, node_id); dp_end[index_i] = GET_DP_END(graph, w, node_id);

        beg = dp_beg[index_i]; end = dp_end[index_i];
        // init h, e
        for (q_i = beg; q_i <= end; ++q_i) {
            dp_h[q_i] = dp_e[q_i] = INT32_MIN;
        } 
        
        for (i = 0; i < pre_n[index_i]; ++i) {
            pre_i = pre_index[index_i][i];
            pre_dp_h = DP_H + pre_i * matrix_col_n; pre_dp_e = DP_E + pre_i * matrix_col_n;
            pre_beg = dp_beg[pre_i]; pre_end = dp_end[pre_i];
            // set M from (pre_i, q_i-1)
            _beg = MAX_OF_TWO(beg-1, pre_beg), _end = MIN_OF_TWO(end-1, pre_end);
            for (q_i = _beg; q_i <= _end; ++q_i) { // SIMD parallelization
                m_pre_i[q_i+1] = pre_dp_h[q_i] > dp_h[q_i+1] ? pre_i : m_pre_i[q_i+1];
                dp_h[q_i+1] = MAX_OF_TWO(pre_dp_h[q_i], dp_h[q_i+1]);
            }
            _beg = MAX_OF_TWO(beg, pre_beg), _end = MIN_OF_TWO(end, pre_end);
            // set E from (pre_i, q_i)
            for (q_i = _beg; q_i <= _end; ++q_i) { // SIMD parallelization
                e_pre_i[q_i] = pre_dp_e[q_i] > dp_e[q_i] ? pre_i : e_pre_i[q_i];
                dp_e[q_i] = MAX_OF_TWO(pre_dp_e[q_i], dp_e[q_i]);
            }
        }
        // compare M, E, and F
        for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
            // get M score
            if (q_i) {
                dp_h[q_i] += q[q_i-1]; mx[q_i-1] = q[q_i-1] == abpt->match ? m0 : x3;
                // h,hd for current cell
                hd[q_i-1] = dp_h[q_i] >= dp_e[q_i] ? (mx[q_i-1] << HOP_OFF_SET | m_pre_i[q_i] << HID_OFF_SET) : (he | e_pre_i[q_i] << HID_OFF_SET);
            }
            dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_e[q_i]);
        }

        // set F from (index_i, q_i-1)
        dp_f[beg] = INT32_MIN;
        for (q_i = beg+1; q_i <= end; ++q_i) { // XXX no SIMD parallelization
            dp_f[q_i] = MAX_OF_TWO(dp_h[q_i-1] - gap_oe, dp_f[q_i-1] - gap_e);
            // fd = h - oe > f - e ? 0 : f2
        }
            
        // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell
        for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
            // h,hd for current cell
            if (q_i) hd[q_i-1] = dp_h[q_i] >= dp_f[q_i] ? hd[q_i-1] : hf;
            dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_f[q_i]);

            // e,ed for next cell
            tmp = dp_h[q_i] - gap_oe;
            dp_e[q_i] -= gap_e;
            if (q_i) ed = dp_e[q_i] > tmp ? (ee | e_pre_i[q_i] << EID_OFF_SET) : (mx[q_i-1] << EOP_OFF_SET | m_pre_i[q_i] << EID_OFF_SET);
            dp_e[q_i] = MAX_OF_TWO(dp_e[q_i], tmp);

            // fd for next cell
            if (q_i) {
                dp_f[q_i] -= gap_e;
                fd = dp_f[q_i] > tmp ? ff : (mx[q_i-1] << FOP_OFF_SET);

                z[q_i-1] = hd[q_i-1] | ed | fd;
            }
        }
    }


#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("index: %d\t", j);
        //for (i = 0; i <= qlen; ++i) {
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {
            printf("%d:(%d,%d)\t", i, DP_H[j*matrix_col_n+i], DP_E[j*matrix_col_n+i]);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    // for global alignment, find best backtrack position
    // TODO semi-global, local
    for (i = 0; i < graph->node[ABPOA_SINK_NODE_ID].in_edge_n; ++i) {
        in_id = graph->node[ABPOA_SINK_NODE_ID].in_id[i];
        in_index = abpoa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, DP_H[in_index * matrix_col_n + qlen], in_index, qlen);
    }

    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
    // backtrack from best score
    abpoa_backtrack(backtrack_z, best_i, best_j, z_col_n, graph, n_cigar, graph_cigar);
    
    free(DP_H); free(DP_E); free(dp_f); free(m_pre_i); free(e_pre_i); free(hd); free(mx);
    free(qp); free(backtrack_z);
    free(dp_beg); free(dp_end);
    for (i = 0; i < graph->node_n; ++i) { 
        free(pre_index[i]); free(next_index[i]); 
    } 
    free(pre_index); free(pre_n); free(next_index); free(next_n);
    return best_score;

    return 0;
}*/

// gap_open1 == 0: ed/fd = hd
int abpoa_banded_global_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar) {
    int **pre_index, *pre_n, pre_i;
    int target_node_n = graph->node_n - 2, matrix_row_n = graph->node_n, matrix_col_n = qlen + 1, z_col_n = qlen;
    int i, j, k, w, *dp_beg, *dp_end, node_id, index_i, q_i;
    int *DP_H, *DP_E, *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *dp_f, tmp; // score type: 8/16/32

    uint8_t *backtrack_z=NULL, *z; // backtrack cell: f<<4|e<<2|h, MATCH:0, DELETION:1, INSERTION:2
    uint8_t *hd=NULL, fd=-1, ed=-1, m0=0x0, e1=0x1, f2=0x2, he, hf, hm, ee, em, ff, fm; 
    int best_score = abpt->inf_min, inf_min = abpt->inf_min, best_i=0, best_j=0;
    int *qp, *mat = abpt->mat, gap_o = abpt->gap_open1, gap_e = abpt->gap_ext1, gap_oe = abpt->gap_open1 + abpt->gap_ext1;

    {   // allocate memory 
        qp = (int*)_err_malloc(qlen * abpt->m * sizeof(int)); // qp should has the same type to H/E/F
        DP_H = (int*)_err_malloc(matrix_row_n * matrix_col_n * sizeof(int));
        DP_E = (int*)_err_malloc(matrix_row_n * matrix_col_n * sizeof(int));
        dp_f = (int*)_err_malloc(matrix_col_n * sizeof(int));

        dp_beg = (int*)_err_malloc((matrix_row_n) * sizeof(int)); dp_end = (int*)_err_calloc(matrix_row_n, sizeof(int));
        // when w <= 0, do whole global
        w = abpt->bw <= 0 ? qlen : abpt->bw;

        if (graph_cigar && n_cigar) {
            backtrack_z = (uint8_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint8_t));

            // for backtrack
            hd = (uint8_t*)_err_malloc(qlen * sizeof(uint8_t));
            hm = m0 << HOP_OFF_SET, he = e1 << HOP_OFF_SET, hf = f2 << HOP_OFF_SET;
            em = m0 << EOP_OFF_SET, ee = e1 << EOP_OFF_SET; fm = m0 << FOP_OFF_SET, ff = f2 << FOP_OFF_SET; 
        }

        pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
        pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));
    }

    {   // generate the query profile
        for (k = i = 0; k < abpt->m; ++k) { // SIMD parallelization
            const int *p = &mat[k * abpt->m];
            for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
        }

        // index of pre-node
        for (i = 0; i < graph->node_n; ++i) {
            node_id = abpoa_graph_index_to_node_id(graph, i); // i: node index
            pre_n[i] = graph->node[node_id].in_edge_n;
            pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
            for (j = 0; j < pre_n[i]; ++j) {
                pre_index[i][j] = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
            }
        }
    }

    {   // DP loop
        // DP cell: H[i,j], E[i+1,j], F[i,j+1]
        // fill the first row
        dp_beg[0] = GET_DP_BEGIN(graph, w, 0, qlen); dp_end[0] = GET_DP_END(graph, w, 0, qlen);

        DP_H[0] = 0; DP_E[0] = -gap_oe;
        for (i = 1; i <= dp_end[0]; ++i) { // SIMD parallelization
            DP_H[i] = -gap_o - gap_e * i; DP_E[i] = inf_min;
        }

        int beg, end, _beg, _end, pre_beg, pre_end;

        for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {
            node_id = abpoa_graph_index_to_node_id(graph, index_i);

            int *q = &qp[graph->node[node_id].base * qlen];
            dp_h = DP_H + index_i * matrix_col_n; dp_e = DP_E + index_i * matrix_col_n;
            z = &backtrack_z[(index_i-1) * qlen];

            dp_beg[index_i] = GET_DP_BEGIN(graph, w, node_id, qlen); dp_end[index_i] = GET_DP_END(graph, w, node_id, qlen);

            beg = dp_beg[index_i]; end = dp_end[index_i];
            // first column
            dp_f[beg] = inf_min;
            if (beg == 0) {
                int min_rank = abpoa_graph_node_id_to_min_rank(graph, node_id);
                dp_h[0] = -(gap_o + gap_e * min_rank);
                dp_e[0] = dp_h[0] - gap_e;
                beg = 1;
            }

            // loop query
            // init h, e                       SIMD parallelization
            for (q_i = beg; q_i <= end; ++q_i) dp_h[q_i] = dp_e[q_i] = inf_min;

            // get max m and e
            for (i = 0; i < pre_n[index_i]; ++i) {
                pre_i = pre_index[index_i][i];
                pre_dp_h = DP_H + pre_i * matrix_col_n; pre_dp_e = DP_E + pre_i * matrix_col_n;
                pre_beg = dp_beg[pre_i]; pre_end = dp_end[pre_i];
                // set M from (pre_i, q_i-1)
                _beg = MAX_OF_TWO(beg-1, pre_beg), _end = MIN_OF_TWO(end-1, pre_end);
                for (q_i = _beg; q_i <= _end; ++q_i) { // SIMD parallelization
                    dp_h[q_i+1] = MAX_OF_TWO(pre_dp_h[q_i], dp_h[q_i+1]);
                }
                // set E from (pre_i, q_i)
                _beg = MAX_OF_TWO(beg, pre_beg), _end = MIN_OF_TWO(end, pre_end);
                for (q_i = _beg; q_i <= _end; ++q_i) { // SIMD parallelization
                    dp_e[q_i] = MAX_OF_TWO(pre_dp_e[q_i], dp_e[q_i]);
                }
            }

            if (n_cigar && graph_cigar) {
                // compare M, E, and F
                for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
                    // get M score; h, hd for current cell
                    dp_h[q_i] += q[q_i-1]; 
                    hd[q_i-1] = dp_h[q_i] >= dp_e[q_i] ? hm : he;
                    dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_e[q_i]);
                }

                // set F from (index_i, q_i-1)
                _beg = dp_beg[index_i] == 0 ? 0 : beg;
                for (q_i = _beg+1; q_i <= end; ++q_i) { // no SIMD parallelization
                    dp_f[q_i] = MAX_OF_TWO(dp_h[q_i-1] - gap_oe, dp_f[q_i-1] - gap_e);
                    // fd = h - oe > f - e ? 0 : f2
                }

                // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell
                for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
                    // h, hd for current cell
                    hd[q_i-1] = dp_h[q_i] >= dp_f[q_i] ? hd[q_i-1] : hf;
                    dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_f[q_i]);

                    // e, ed for next cell
                    tmp = dp_h[q_i] - gap_oe;
                    dp_e[q_i] -= gap_e;
                    ed = gap_o ? (dp_e[q_i] > tmp ? ee : em) : (hd[q_i-1] == hm ? em : ee);
                    dp_e[q_i] = MAX_OF_TWO(dp_e[q_i], tmp);

                    // fd for next cell
                    dp_f[q_i] -= gap_e;
                    fd = gap_o ? (dp_f[q_i] > tmp ? ff : fm) : (hd[q_i-1] == hm ? fm : ff);
                    z[q_i-1] = hd[q_i-1] | ed | fd;
                }
            } else {
                // compare M, E, and F
                for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
                    // get M score; h, hd for current cell
                    dp_h[q_i] += q[q_i-1]; 
                    dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_e[q_i]);
                }

                // set F from (index_i, q_i-1)
                _beg = dp_beg[index_i] == 0 ? 0 : beg;
                for (q_i = _beg+1; q_i <= end; ++q_i) { // no SIMD parallelization
                    dp_f[q_i] = MAX_OF_TWO(dp_h[q_i-1] - gap_oe, dp_f[q_i-1] - gap_e);
                    // fd = h - oe > f - e ? 0 : f2
                }

                // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell
                for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
                    // h, hd for current cell
                    dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_f[q_i]);

                    // e, ed for next cell
                    tmp = dp_h[q_i] - gap_oe;
                    dp_e[q_i] -= gap_e;
                    dp_e[q_i] = MAX_OF_TWO(dp_e[q_i], tmp);

                    // fd for next cell
                    dp_f[q_i] -= gap_e;
                    z[q_i-1] = hd[q_i-1] | ed | fd;
                }
            }
        }
    }

#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("index: %d\t", j);
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {
            printf("%d:(%d,%d)\t", i, DP_H[j*matrix_col_n+i], DP_E[j*matrix_col_n+i]);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    for (i = 0; i < graph->node[ABPOA_SINK_NODE_ID].in_edge_n; ++i) { // for global alignment, find best backtrack position
        in_id = graph->node[ABPOA_SINK_NODE_ID].in_id[i];
        in_index = abpoa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, DP_H[in_index * matrix_col_n + qlen], in_index, qlen);
    }

#ifdef __DEBUG__
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
#endif
    { // backtrack from best score
        if (n_cigar && graph_cigar) abpoa_backtrack(DP_H, DP_E, matrix_col_n, abpt->m, abpt->mat, abpt->gap_ext1, pre_index, pre_n, backtrack_z, best_i, best_j, z_col_n, graph, query, n_cigar, graph_cigar);
    }

    { // free variables
        free(DP_H); free(DP_E); free(dp_f); free(qp); free(dp_beg); free(dp_end);
        if (n_cigar && graph_cigar) free(backtrack_z), free(hd);
        for (i = 0; i < graph->node_n; ++i) free(pre_index[i]); 
        free(pre_index); free(pre_n);
    }
    return best_score;
}

int ada_abpoa_banded_global_align_sequence_with_graph(abpoa_graph_t *graph, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar) {
    int **pre_index, *pre_n, pre_i;
    int target_node_n = graph->node_n - 2, matrix_row_n = graph->node_n, matrix_col_n = qlen + 1, z_col_n = qlen;
    int i, j, k, w, *dp_beg, *dp_end, node_id, index_i, q_i;
    int *DP_H, *DP_E, *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *dp_f, tmp; // score type: 8/16/32

    uint8_t *backtrack_z=NULL, *z; // backtrack cell: f<<4|e<<2|h, MATCH:0, DELETION:1, INSERTION:2
    uint8_t *hd=NULL, fd=-1, ed=-1, m0=0x0, e1=0x1, f2=0x2, he, hf, hm, ee, em, ff, fm; 
    int best_score = abpt->inf_min, inf_min = abpt->inf_min, best_i=0, best_j=0;
    int *qp, *mat = abpt->mat, gap_o = abpt->gap_open1, gap_e = abpt->gap_ext1, gap_oe = abpt->gap_open1 + abpt->gap_ext1;

    {   // allocate memory 
        qp = (int*)_err_malloc(qlen * abpt->m * sizeof(int)); // qp should has the same type to H/E/F
        DP_H = (int*)_err_malloc(matrix_row_n * matrix_col_n * sizeof(int));
        DP_E = (int*)_err_malloc(matrix_row_n * matrix_col_n * sizeof(int));
        dp_f = (int*)_err_malloc(matrix_col_n * sizeof(int));

        dp_beg = (int*)_err_malloc((matrix_row_n) * sizeof(int)); dp_end = (int*)_err_calloc(matrix_row_n, sizeof(int));
        // when w <= 0, do whole global
        w = abpt->bw <= 0 ? qlen : abpt->bw;

        if (graph_cigar && n_cigar) {
            backtrack_z = (uint8_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint8_t));

            // for backtrack
            hd = (uint8_t*)_err_malloc(qlen * sizeof(uint8_t));
            hm = m0 << HOP_OFF_SET, he = e1 << HOP_OFF_SET, hf = f2 << HOP_OFF_SET;
            em = m0 << EOP_OFF_SET, ee = e1 << EOP_OFF_SET; fm = m0 << FOP_OFF_SET, ff = f2 << FOP_OFF_SET; 
        }

        pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
        pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));
    }

    {   // generate the query profile
        for (k = i = 0; k < abpt->m; ++k) { // SIMD parallelization
            const int *p = &mat[k * abpt->m];
            for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
        }

        // index of pre-node
        for (i = 0; i < graph->node_n; ++i) {
            node_id = abpoa_graph_index_to_node_id(graph, i); // i: node index
            pre_n[i] = graph->node[node_id].in_edge_n;
            pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
            for (j = 0; j < pre_n[i]; ++j) {
                pre_index[i][j] = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
            }
        }

        // init min/max_rank
        for (i = 2; i < graph->node_n; ++i) {
            node_id = graph->index_to_node_id[i];
            graph->node_id_to_min_rank[node_id] = graph->node_n;
            graph->node_id_to_max_rank[node_id] = 0;
        }
    }

    {   // DP loop
        // DP cell: H[i,j], E[i+1,j], F[i,j+1]
        int beg, end, _beg, _end, pre_beg, pre_end;
        // fill the first row
        // get beg and end based on min/max_rank and min/max_remain
        dp_beg[0] = GET_DP_BEGIN(graph, w, 0, qlen); dp_end[0] = GET_DP_END(graph, w, 0, qlen);

        DP_H[0] = 0; DP_E[0] = -gap_oe;
        for (i = 1; i <= dp_end[0]; ++i) { // SIMD parallelization
            DP_H[i] = -gap_o - gap_e * i; DP_E[i] = inf_min;
        }

        for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {
            node_id = abpoa_graph_index_to_node_id(graph, index_i);

            int *q = &qp[graph->node[node_id].base * qlen];
            dp_h = DP_H + index_i * matrix_col_n; dp_e = DP_E + index_i * matrix_col_n;
            z = &backtrack_z[(index_i-1) * qlen];

            //printf("index: %d, min_rank: %d, max_rank: %d\n", index_i, graph->node_id_to_min_rank[node_id], graph->node_id_to_max_rank[node_id]);
            dp_beg[index_i] = GET_DP_BEGIN(graph, w, node_id, qlen); dp_end[index_i] = GET_DP_END(graph, w, node_id, qlen);

            beg = dp_beg[index_i]; end = dp_end[index_i];
            // first column
            dp_f[beg] = inf_min;
            if (beg == 0) {
                int min_rank = abpoa_graph_node_id_to_min_rank(graph, node_id);
                dp_h[0] = -(gap_o + gap_e * min_rank);
                dp_e[0] = dp_h[0] - gap_e;
                beg = 1;
            }

            // loop query
            // init h, e                       SIMD parallelization
            for (q_i = beg; q_i <= end; ++q_i) dp_h[q_i] = dp_e[q_i] = inf_min;

            // get max m and e
            for (i = 0; i < pre_n[index_i]; ++i) {
                pre_i = pre_index[index_i][i];
                pre_dp_h = DP_H + pre_i * matrix_col_n; pre_dp_e = DP_E + pre_i * matrix_col_n;
                pre_beg = dp_beg[pre_i]; pre_end = dp_end[pre_i];
                // set M from (pre_i, q_i-1)
                _beg = MAX_OF_TWO(beg-1, pre_beg), _end = MIN_OF_TWO(end-1, pre_end);
                for (q_i = _beg; q_i <= _end; ++q_i) { // SIMD parallelization
                    dp_h[q_i+1] = MAX_OF_TWO(pre_dp_h[q_i], dp_h[q_i+1]);
                }
                // set E from (pre_i, q_i)
                _beg = MAX_OF_TWO(beg, pre_beg), _end = MIN_OF_TWO(end, pre_end);
                for (q_i = _beg; q_i <= _end; ++q_i) { // SIMD parallelization
                    dp_e[q_i] = MAX_OF_TWO(pre_dp_e[q_i], dp_e[q_i]);
                }
            }

            if (n_cigar && graph_cigar) {
                // compare M, E, and F
                for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
                    // get M score; h, hd for current cell
                    dp_h[q_i] += q[q_i-1]; 
                    hd[q_i-1] = dp_h[q_i] >= dp_e[q_i] ? hm : he;
                    dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_e[q_i]);
                }

                // set F from (index_i, q_i-1)
                _beg = dp_beg[index_i] == 0 ? 0 : beg;
                for (q_i = _beg+1; q_i <= end; ++q_i) { // no SIMD parallelization
                    dp_f[q_i] = MAX_OF_TWO(dp_h[q_i-1] - gap_oe, dp_f[q_i-1] - gap_e);
                    // fd = h - oe > f - e ? 0 : f2
                }

                // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell
                for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
                    // h, hd for current cell
                    hd[q_i-1] = dp_h[q_i] >= dp_f[q_i] ? hd[q_i-1] : hf;
                    dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_f[q_i]);

                    // e, ed for next cell
                    tmp = dp_h[q_i] - gap_oe;
                    dp_e[q_i] -= gap_e;
                    ed = gap_o ? (dp_e[q_i] > tmp ? ee : em) : (hd[q_i-1] == hm ? em : ee);
                    dp_e[q_i] = MAX_OF_TWO(dp_e[q_i], tmp);

                    // fd for next cell
                    dp_f[q_i] -= gap_e;
                    fd = gap_o ? (dp_f[q_i] > tmp ? ff : fm) : (hd[q_i-1] == hm ? fm : ff);
                    z[q_i-1] = hd[q_i-1] | ed | fd;
                }
            } else {
                // compare M, E, and F
                for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
                    // get M score; h, hd for current cell
                    dp_h[q_i] += q[q_i-1]; 
                    dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_e[q_i]);
                }

                // set F from (index_i, q_i-1)
                _beg = dp_beg[index_i] == 0 ? 0 : beg;
                for (q_i = _beg+1; q_i <= end; ++q_i) { // no SIMD parallelization
                    dp_f[q_i] = MAX_OF_TWO(dp_h[q_i-1] - gap_oe, dp_f[q_i-1] - gap_e);
                    // fd = h - oe > f - e ? 0 : f2
                }

                // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell
                for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
                    // h, hd for current cell
                    dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_f[q_i]);

                    // e, ed for next cell
                    tmp = dp_h[q_i] - gap_oe;
                    dp_e[q_i] -= gap_e;
                    dp_e[q_i] = MAX_OF_TWO(dp_e[q_i], tmp);

                    // fd for next cell
                    dp_f[q_i] -= gap_e;
                    z[q_i-1] = hd[q_i-1] | ed | fd;
                }
            }

            // XXX need hd here
            { // set current rank, set min/max rank for next nodes
                // select max dp_h
                int max = INT32_MIN, max_i = -1;
                for (q_i = beg; q_i <= end; ++q_i) { // beg >= 1
                    if (dp_h[q_i] > max) {
                        max_i = q_i;
                        max = dp_h[q_i];
                    }
                }
                // determine max pre_i, then determin current rank
                int which, s, max_pre_i=-1, max_pre_id;
                which = (hd[max_i-1] >> HOP_OFF_SET) & 3;
                if (which == 0) { // match
                    s = graph->node[node_id].base == query[max_i-1] ? abpt->match : -abpt->mismatch;
                    for (k = 0; k < pre_n[index_i]; ++k) {
                        pre_i = pre_index[index_i][k];
                        if (DP_H[pre_i * matrix_col_n + max_i-1] + s == DP_H[index_i * matrix_col_n + max_i]) {
                            max_pre_i = pre_i;
                            break;
                        }
                    }
                    max_pre_id = abpoa_graph_index_to_node_id(graph, max_pre_i);
                } else if (which == 1) { // deletion
                    for (k = 0; k < pre_n[index_i]; ++k) {
                        pre_i = pre_index[index_i][k];
                        if (DP_E[pre_i * matrix_col_n + max_i] == DP_H[index_i * matrix_col_n + max_i]) {
                            max_pre_i = pre_i;
                            break;
                        }
                    }
                    max_pre_id = abpoa_graph_index_to_node_id(graph, max_pre_i);
                } else { // insertion
                    err_fatal_simple("Unexpected cigar op.");
                }  
                // set min/max_rank for next nodes
                graph->node_id_to_max_rank[node_id] = graph->node_id_to_max_rank[max_pre_id] + 1;
                graph->node_id_to_min_rank[node_id] = graph->node_id_to_min_rank[max_pre_id] + 1;
                int out_rank = graph->node_id_to_max_rank[node_id] + 1;
                for (i = 0; i < graph->node[node_id].out_edge_n; ++i) {
                    int out_node_id = graph->node[node_id].out_id[i];
                    if (out_rank > graph->node_id_to_max_rank[out_node_id])
                        graph->node_id_to_max_rank[out_node_id] = out_rank;
                    if (out_rank < graph->node_id_to_min_rank[out_node_id])
                        graph->node_id_to_min_rank[out_node_id] = out_rank;
                }
            }
        }
    }

#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("index: %d\t", j);
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {
            printf("%d:(%d,%d)\t", i, DP_H[j*matrix_col_n+i], DP_E[j*matrix_col_n+i]);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    for (i = 0; i < graph->node[ABPOA_SINK_NODE_ID].in_edge_n; ++i) { // for global alignment, find best backtrack position
        in_id = graph->node[ABPOA_SINK_NODE_ID].in_id[i];
        in_index = abpoa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, DP_H[in_index * matrix_col_n + qlen], in_index, qlen);
    }

#ifdef __DEBUG__
    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
#endif
    { // backtrack from best score
        if (n_cigar && graph_cigar) abpoa_backtrack(DP_H, DP_E, matrix_col_n, abpt->m, abpt->mat, abpt->gap_ext1, pre_index, pre_n, backtrack_z, best_i, best_j, z_col_n, graph, query, n_cigar, graph_cigar);
    }

    { // free variables
        free(DP_H); free(DP_E); free(dp_f); free(qp); free(dp_beg); free(dp_end);
        if (n_cigar && graph_cigar) free(backtrack_z), free(hd);
        for (i = 0; i < graph->node_n; ++i) free(pre_index[i]); 
        free(pre_index); free(pre_n);
    }
    return best_score;
}
// TODO local, extend, semi-global
