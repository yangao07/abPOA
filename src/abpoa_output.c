#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "abpoa_output.h"
#include "abpoa.h"
#include "abpoa_graph.h"
#include "utils.h"
#include "abpoa_seq.h"
#include "kdq.h"

extern char ab_char256_table[256];
char ab_LogTable65536[65536];
char ab_bit_table16[65536];

#define NAT_E 2.718281828459045
static const char ab_LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
    uint32_t t, tt;
    if ((tt = v>>16)) return (t = tt>>8) ? 24 + ab_LogTable256[t] : 16 + ab_LogTable256[tt];
    return (t = v>>8) ? 8 + ab_LogTable256[t] : ab_LogTable256[v];
}

void set_65536_table(void) {
    int i;
    for (i = 0; i < 65536; ++i) {
        ab_LogTable65536[i] = ilog2_32(i);
    }
}

void set_bit_table16(void) {
    int i; ab_bit_table16[0] = 0;
    for (i = 0; i != 65536; ++i) ab_bit_table16[i] = (i&1) + ab_bit_table16[i>>1];
}

#define get_bit_cnt4(table, b) (table[(b)&0xffff] + table[(b)>>16&0xffff] + table[(b)>>32&0xffff] + table[(b)>>48&0xffff])

static inline int ilog2_64(uint64_t v) {
    uint64_t t, tt;
    if ((tt = v >> 32)) return (t = tt >> 16) ? 48 + ab_LogTable65536[t] : 32 + ab_LogTable65536[tt];
    return (t = v>>16) ? 16 + ab_LogTable65536[t] : ab_LogTable65536[v];
}

KDQ_INIT(int)
#define kdq_int_t kdq_t(int)

static inline int get_read_cnt(uint64_t *read_ids, int read_ids_n) {
    int i, c;
    for (i = c =0; i < read_ids_n; ++i) {
        c += get_bit_cnt4(ab_bit_table16, read_ids[i]);
    }
    return c;
}

abpoa_cons_t *abpoa_allocate_rc_msa(abpoa_cons_t *abc, int msa_len, int n_seq, int n_cons) {
    int i;
    abc->n_seq = n_seq; abc->msa_len = msa_len;
    abc->msa_base = (uint8_t**)_err_malloc((n_seq+n_cons) * sizeof(uint8_t*));
    for (i = 0; i < n_seq+n_cons; ++i) {
        abc->msa_base[i] = (uint8_t*)_err_malloc(msa_len * sizeof(uint8_t));
    }
    return abc;
}

void abpoa_output_rc_msa(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp) {
    if (out_fp == NULL) return;
    int i, j;
    abpoa_seq_t *abs = ab->abs; abpoa_cons_t *abc = ab->abc;
    if (abc->msa_len <= 0) return;
    for (i = 0; i < abs->n_seq; ++i) {
        if (abs->name[i].l > 0) {
            if (abs->is_rc[i]) fprintf(out_fp, ">%s_reverse_complement\n", abs->name[i].s);
            else fprintf(out_fp, ">%s\n", abs->name[i].s);
        } else {
            fprintf(out_fp, ">Seq_%d\n", i+1);
        }
        for (j = 0; j < abc->msa_len; ++j) fprintf(out_fp, "%c", ab_char256_table[abc->msa_base[i][j]]);
        fprintf(out_fp, "\n");
    }
    if (abpt->out_cons) { // RC-MSA for consensus sequence
        int cons_i;
        for (cons_i = 0; cons_i < abc->n_cons; cons_i++) {
            fprintf(out_fp, ">Consensus_sequence");
            if (abc->n_cons > 1) {
                fprintf(out_fp, "_%d ", cons_i+1);
                for (j = 0; j < abc->clu_n_seq[cons_i]; ++j) { // cluter read_id
                    if (j != 0) fprintf(out_fp, ",");
                    fprintf(out_fp, "%d", abc->clu_read_ids[cons_i][j]);
                }
            }
            fprintf(out_fp, "\n");
            for (i = 0; i < abc->msa_len; ++i) fprintf(out_fp, "%c", ab_char256_table[abc->msa_base[abc->n_seq+cons_i][i]]);
            fprintf(out_fp, "\n");
        }
    }
}

void abpoa_set_msa_seq(abpoa_node_t node, int rank, uint8_t **msa_base) {
    int i, j, b, read_id; uint8_t base = node.base;
    uint64_t num, tmp;

    b = 0;
    for (i = 0; i < node.read_ids_n; ++i) {
        for (j = 0; j < node.out_edge_n; ++j) {
            num = node.read_ids[j][i];
            while (num) {
                tmp = num & -num;
                read_id = ilog2_64(tmp);
                msa_base[b+read_id][rank-1] = base;
                num ^= tmp;
            }
        }
        b += 64;
    }
}

int abpoa_collect_msa(abpoa_graph_t *abg, abpoa_para_t *abpt, uint8_t **msa, int n_seq) {
    if (abg->node_n <= 2) return 0;
    abpoa_set_msa_rank(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    int msa_len = abg->node_id_to_msa_rank[ABPOA_SINK_NODE_ID]-1;
    int i, j;
    for (i = 0; i < n_seq; ++i) {
        msa[i] = (uint8_t*)_err_malloc(msa_len * sizeof(uint8_t));
        for (j = 0; j < msa_len; ++j) msa[i][j] = abpt->m;
    }
    int rank, aligned_id;
    // if (out_fp && abpt->out_msa_header == 0) fprintf(out_fp, ">Multiple_sequence_alignment\n");
    for (i = 2; i < abg->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(abg, i);
        for (j = 0; j < abg->node[i].aligned_node_n; ++j) {
            aligned_id = abg->node[i].aligned_node_id[j];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
        }
        // assign seq
        abpoa_set_msa_seq(abg->node[i], rank, msa);
    }
    return msa_len;
}

// only generate rc-msa, output in separated func
void abpoa_generate_rc_msa(abpoa_t *ab, abpoa_para_t *abpt) {
    abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) return;
    abpoa_set_msa_rank(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    if (abpt->out_cons) abpoa_generate_consensus(ab, abpt);

    abpoa_seq_t *abs = ab->abs; abpoa_cons_t *abc = ab->abc;
    int i, j, aligned_id, n_seq = abs->n_seq;
    int msa_len = abg->node_id_to_msa_rank[ABPOA_SINK_NODE_ID]-1;

    abpoa_allocate_rc_msa(abc, msa_len, n_seq, abc->n_cons);
    for (i = 0; i < n_seq; ++i) {
        for (j = 0; j < abc->msa_len; ++j) 
            abc->msa_base[i][j] = abpt->m;
    }

    int rank;
    // if (out_fp && abpt->out_msa_header == 0) fprintf(out_fp, ">Multiple_sequence_alignment\n");
    for (i = 2; i < abg->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(abg, i);
        for (j = 0; j < abg->node[i].aligned_node_n; ++j) {
            aligned_id = abg->node[i].aligned_node_id[j];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
        }
        // assign seq
        abpoa_set_msa_seq(abg->node[i], rank, abc->msa_base);
    }
    if (abpt->out_cons) {
        int cons_i, cur_id;
        for (cons_i = 0; cons_i < abc->n_cons; cons_i++) {
            for (i = 0; i < msa_len; ++i) abc->msa_base[n_seq+cons_i][i] = abpt->m;
            for (i = 0; i < abc->cons_len[cons_i]; ++i) {
                cur_id = abc->cons_node_ids[cons_i][i];
                rank = abpoa_graph_node_id_to_msa_rank(abg, cur_id);
                for (j = 0; j < abg->node[cur_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[cur_id].aligned_node_id[j];
                    rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
                }
                abc->msa_base[n_seq+cons_i][rank-1] = abc->cons_base[cons_i][i];
            }
        }
    }
}

// generate & output gfa
void abpoa_generate_gfa(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp) {
    if (out_fp == NULL) return;
    abpoa_seq_t *abs = ab->abs; abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) return;

    // traverse graph 
    int *in_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    int n_seq = abs->n_seq;
    int **read_paths = (int**)_err_malloc(n_seq * sizeof(int*)), *read_path_i = (int*)_err_calloc(n_seq, sizeof(int));
    int i, j, cur_id, pre_id, out_id, *id;
    for (i = 0; i < abg->node_n; ++i) in_degree[i] = abg->node[i].in_edge_n;
    for (i = 0; i < n_seq; ++i) read_paths[i] = (int*)_err_malloc(abg->node_n * sizeof(int));

    // output comment and header
    int nl = 0;
    for (i = 2; i < abg->node_n; ++i) nl += abg->node[i].in_edge_n;
    fprintf(out_fp, "H\tVN:Z:1.0\tNS:i:%d\tNL:i:%d\tNP:i:%d\n", abg->node_n-2, nl - abg->node[ABPOA_SRC_NODE_ID].out_edge_n, n_seq + abpt->out_cons);

    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, ABPOA_SRC_NODE_ID); 
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == ABPOA_SINK_NODE_ID) {
            kdq_destroy_int(q);
            break;
        } else {
            if (cur_id != ABPOA_SRC_NODE_ID) {
                // output node
                fprintf(out_fp, "S\t%d\t%c\n", cur_id-1, ab_char256_table[abg->node[cur_id].base]);
                // output all links based pre_ids
                for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
                    pre_id = abg->node[cur_id].in_id[i];
                    if (pre_id != ABPOA_SRC_NODE_ID)
                        fprintf(out_fp, "L\t%d\t+\t%d\t+\t0M\n", pre_id-1, cur_id-1);
                }
                // add node id to read path
                int b, read_id; uint64_t num, tmp;
                b = 0;
                for (i = 0; i < abg->node[cur_id].read_ids_n; ++i) {
                    for (j = 0; j < abg->node[cur_id].out_edge_n; ++j) {
                        num = abg->node[cur_id].read_ids[j][i];
                        while (num) {
                            tmp = num & -num;
                            read_id = ilog2_64(tmp);
                            read_paths[b+read_id][read_path_i[b+read_id]++] = cur_id-1;
                            num ^= tmp;
                        }
                    }
                    b += 64;
                }
            }
            for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                out_id = abg->node[cur_id].out_id[i];
                if (--in_degree[out_id] == 0) {
                    kdq_push_int(q, out_id);
                }
            }
        }
    }
    // output read paths
    for (i = 0; i < n_seq; ++i) {
        if (abs->name[i].l > 0) fprintf(out_fp, "P\t%s\t", abs->name[i].s);
        else fprintf(out_fp, "P\t%d\t", i+1);
        if (abs->is_rc[i]) {
            for (j = read_path_i[i]-1; j >= 0; --j) {
                fprintf(out_fp, "%d-", read_paths[i][j]);
                if (j != 0) fprintf(out_fp, ",");
                else fprintf(out_fp, "\t*\n");
            }
        } else {
            for (j = 0; j < read_path_i[i]; ++j) {
                fprintf(out_fp, "%d+", read_paths[i][j]);
                if (j != read_path_i[i]-1) fprintf(out_fp, ",");
                else fprintf(out_fp, "\t*\n");
            }
        }
    }
    if (abpt->out_cons) {
        abpoa_generate_consensus(ab, abpt);
        abpoa_cons_t *abc = ab->abc;
        int cons_i;
        for (cons_i = 0; cons_i < abc->n_cons; ++cons_i) {
            fprintf(out_fp, "P\tConsensus_sequence");
            if (abc->n_cons > 1) fprintf(out_fp, "_%d", cons_i+1);
            fprintf(out_fp, "\t");
            for (i = 0; i < abc->cons_len[cons_i]; ++i) {
                cur_id = abc->cons_node_ids[cons_i][i];
                fprintf(out_fp, "%d+", cur_id-1);
                if (i != abc->cons_len[cons_i]-1) fprintf(out_fp, ",");
                else fprintf(out_fp, "\t*\n"); 
            }

        }
    }
    free(in_degree);
    for (i = 0; i < n_seq; ++i) free(read_paths[i]); 
    free(read_paths); free(read_path_i);
}

int abpoa_cons_phred_score(int n_cov, int n_seq) {
    if (n_cov > n_seq) err_fatal(__func__, "Error: unexpected n_cov/n_seq (%d/%d).", n_cov, n_seq);
    double x, p;
    x = 13.8 * (1.25 * n_cov / n_seq - 0.25);
    p = 1 - 1.0 / (1.0 + pow(NAT_E, -1 * x));
    return (33 + (int)(-10 * log10(p) + 0.499));
}

int get_edge_inclu_read_count(int edge_i, int cons_i, abpoa_node_t *node, uint64_t **clu_read_ids) {
    int n = 0, i; uint64_t b;
    for (i = 0; i < node->read_ids_n; ++i) {
        b = node->read_ids[edge_i][i] & clu_read_ids[cons_i][i];
        n += get_bit_cnt4(ab_bit_table16, b);
    }
    return n;
}

// for n_clu >= 2, get weight of reads within cons_i'th clu
int get_inclu_edge_weight(int edge_i, int cons_i, abpoa_node_t *node, uint64_t **clu_read_ids, int use_qv) {
    // collect read-wise weight
    if (use_qv == 0) return get_edge_inclu_read_count(edge_i, cons_i, node, clu_read_ids);
    int w = 0, i; uint64_t c, one = 1;
    for (i = 0; i < node->m_read; ++i) {
        if (node->read_weight[i] > 0) {
            int n = i / 64, b = i & 0x3f;
            c = node->read_ids[edge_i][n] & clu_read_ids[cons_i][n];
            if ((c & (one << b)) > 0)
                w += node->read_weight[i];
        }
    }
    return w;
}

int get_edge_weight(int edge_i, int cons_i, abpoa_node_t *node, uint64_t **clu_read_ids, int use_qv, int n_clu) {
    if (n_clu == 1) { // use weight directly
        return node->out_edge_weight[edge_i];
    } else { // collect weight of reads within cons_i'th clu
        return get_inclu_edge_weight(edge_i, cons_i, node, clu_read_ids, use_qv);
    }
}

int get_node_weight(int n_clu, int cons_i, abpoa_node_t *node, uint64_t **clu_read_ids, int use_qv) {
    int w = 0, i;
    for (i = 0; i < node->out_edge_n; ++i) {
        w += get_edge_weight(i, cons_i, node, clu_read_ids, use_qv, n_clu);
    }
    return w;
}

// get base coverage for node[id]
int abpoa_node_out_cov(abpoa_node_t *node, int id, uint64_t **clu_read_ids, int cons_i, int n_cons) {
    if (n_cons == 1) return node[id].n_read;
    int i, out_cov = 0;
    for (i = 0; i < node[id].out_edge_n; ++i) {
        out_cov += get_edge_inclu_read_count(i, cons_i, node+id, clu_read_ids);
    }
    return out_cov;
}

int abpoa_node_in_cov(abpoa_node_t *node, int id, uint64_t **clu_read_ids, int cons_i, int n_cons) {
    if (n_cons == 1) return node[id].n_read;
    int i, j, in_id, in_cov = 0;
    // for each id: get max{left_weigth, right_weight}
    for (i = 0; i < node->in_edge_n; ++i) {
        in_id = node[id].in_id[i];
        for (j = 0; j < node[in_id].out_edge_n; ++j) {
            if (node[in_id].out_id[j] == id) {
                in_cov += get_edge_inclu_read_count(j, cons_i, node+in_id, clu_read_ids);
                break;
            }
        }
    }
    return in_cov;
}
int abpoa_node_cov(abpoa_node_t *node, int id, uint64_t **clu_read_ids, int cons_i, int n_cons) {
    if (n_cons == 1) return node[id].n_read;
    return MAX_OF_TWO(abpoa_node_in_cov(node, id, clu_read_ids, cons_i, n_cons), abpoa_node_out_cov(node, id, clu_read_ids, cons_i, n_cons));
}

void abpoa_set_hb_cons(abpoa_graph_t *abg, int **max_out_id, int n_cons, uint64_t **clu_read_ids, int src_id, int sink_id, abpoa_cons_t *abc) {
    abc->n_cons = n_cons;
    int i, j, cur_id;
    for (i = 0; i < n_cons; ++i) {
        cur_id = max_out_id[i][src_id];
        j = 0;
        while (cur_id != sink_id) {
            abc->cons_node_ids[i][j] = cur_id;
            abc->cons_base[i][j] = abg->node[cur_id].base;
            abc->cons_cov[i][j] = abpoa_node_cov(abg->node, cur_id, clu_read_ids, i, n_cons);
            abc->cons_phred_score[i][j] = abpoa_cons_phred_score(abc->cons_cov[i][j], abc->clu_n_seq[i]);
            ++j;
            cur_id = max_out_id[i][cur_id];
        }
        abc->cons_len[i] = j;
    }
}

void abpoa_set_major_voting_cons(abpoa_graph_t *abg, int use_span_read_count, int m, int ***row_column_count, int **msa_node_id, int src_id, int sink_id, int msa_l, abpoa_cons_t *abc) {
    int cur_id, cons_i, i, j, max_c, max_base, total_c, gap_c, c;
    int cons_l;
    for (cons_i = 0; cons_i < abc->n_cons; ++cons_i) {
        cons_l = 0;
        for (i = 0; i < msa_l; ++i) {
            max_c = 0, total_c = 0, max_base = m; //, gap_c = abc->clu_n_seq[cons_i];
            for (j = 0; j < m-1; ++j) {
                c = row_column_count[cons_i][i][j];
                if (c > max_c) {
                    max_c = c;
                    max_base = j;
                }
                total_c += c;
                // gap_c -= c;
            }
            // assert(max_base != m);
            if (use_span_read_count) gap_c = abg->node[msa_node_id[i][max_base]].n_span_read - total_c;
            else gap_c = abc->clu_n_seq[cons_i] - total_c;
            if (max_c >= gap_c) { // append consensus base to abc
                cur_id = msa_node_id[i][max_base];
                abc->cons_node_ids[cons_i][cons_l] = cur_id;
                abc->cons_base[cons_i][cons_l] = max_base;
                abc->cons_cov[cons_i][cons_l] = max_c;
                abc->cons_phred_score[cons_i][cons_l] = abpoa_cons_phred_score(abc->cons_cov[cons_i][cons_l], abc->clu_n_seq[cons_i]);
                cons_l++;
            }
            // fprintf(stderr, "i: %d, max_base: %d, max_c: %d, gap_c: %d, total_c: %d, n_span: %d\n", i, max_base, max_c, gap_c, total_c, abg->node[msa_node_id[i][max_base]].n_span_read);
        }
        abc->cons_len[cons_i] = cons_l;
    }
}

void abpoa_set_row_column_weight(abpoa_graph_t *abg, int n_clu, int m, int ***rc_weight, uint64_t **clu_read_ids, int **msa_node_id) {
    int i, cons_i, k, rank, aligned_id;
    int node_w;
    for (i = 2; i < abg->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(abg, i);
        for (k = 0; k < abg->node[i].aligned_node_n; ++k) {
            aligned_id = abg->node[i].aligned_node_id[k];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
        }
        msa_node_id[rank-1][abg->node[i].base] = i;
        // assign seq
        for (cons_i = 0; cons_i < n_clu; ++cons_i) {
            node_w = abpoa_node_out_cov(abg->node, i, clu_read_ids, cons_i, n_clu);
            rc_weight[cons_i][rank-1][abg->node[i].base] = node_w;
            // for (k = 0; k < abg->node[i].read_ids_n; ++k) {
            //     for (j = 0; j < abg->node[i].out_edge_n; ++j) {
            //         if (n_clu > 1) b = abg->node[i].read_ids[j][k] & clu_read_ids[cons_i][k];
            //         else b = abg->node[i].read_ids[j][k];
            //         rc_count[cons_i][rank-1][abg->node[i].base] += get_bit_cnt4(ab_bit_table16, b);
            //     }
            // }
            rc_weight[cons_i][rank-1][m-1] -= rc_weight[cons_i][rank-1][abg->node[i].base];
        }
    }
}

void set_clu_read_ids(abpoa_cons_t *abc, uint64_t **read_ids, int cons_i, int n_seq) {
    int n, i, j; uint64_t b, one = 1;
    for (i = n = 0; i < n_seq; ++i) {
        j = i / 64; b = i & 0x3f;
        if (read_ids[cons_i][j] & (one << b)) {
            abc->clu_read_ids[cons_i][n++] = i;
        }
    }
    if (n != abc->clu_n_seq[cons_i])
        err_fatal(__func__, "Error in set cluster read ids. (%d, %d)", n, abc->clu_n_seq[cons_i]);
}

// heaviest_bundling
// 1. argmax{cur->weight}
// 2. argmax{out_node->weight}
void abpoa_heaviest_bundling(abpoa_graph_t *abg, abpoa_para_t *abpt, int src_id, int sink_id, int *out_degree, int n_clu, int read_ids_n, uint64_t **clu_read_ids, abpoa_cons_t *abc) {
    int *id, cons_i, i, cur_id, in_id, out_id, max_id; int max_w, out_w;

    int *_out_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    int *score = (int*)_err_malloc(abg->node_n * sizeof(int));
    int **max_out_id = (int**)_err_malloc(n_clu * sizeof(int*));
    for (i = 0; i < n_clu; ++i) max_out_id[i] = (int*)_err_malloc(abg->node_n * sizeof(int));
    if (n_clu == 1) abc->clu_n_seq[0] = abc->n_seq;
    else {
        for (cons_i = 0; cons_i < n_clu; cons_i++) {
            abc->clu_n_seq[cons_i] = get_read_cnt(clu_read_ids[cons_i], read_ids_n);
            set_clu_read_ids(abc, clu_read_ids, cons_i, abc->n_seq);
        }
    }

    for (cons_i = 0; cons_i < n_clu; cons_i++) {
        for (i = 0; i < abg->node_n; ++i) _out_degree[i] = out_degree[i];
        kdq_int_t *q = kdq_init_int();
        kdq_push_int(q, sink_id);
        // reverse Breadth-First-Search
        while ((id = kdq_shift_int(q)) != 0) {
            cur_id = *id;
            if (cur_id == sink_id) {
                max_out_id[cons_i][cur_id] = -1;
                score[cur_id] = 0;
            } else {
                max_id = -1;
                if (cur_id == src_id) {
                    int path_score = -1, path_max_w = -1;
                    for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                        out_id = abg->node[cur_id].out_id[i];
                        out_w = get_edge_weight(i, cons_i, abg->node+cur_id, clu_read_ids, abpt->use_qv, n_clu);
                        if (out_w > path_max_w || (out_w == path_max_w && score[out_id] > path_score)) {
                            max_id = out_id;
                            path_score = score[out_id];
                            path_max_w = out_w;
                        }
                    }
                    max_out_id[cons_i][cur_id] = max_id;
                    kdq_destroy_int(q);
                    break;
                } else {
                    max_w = INT32_MIN;
                    for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                        out_id = abg->node[cur_id].out_id[i];
                        out_w = get_edge_weight(i, cons_i, abg->node+cur_id, clu_read_ids, abpt->use_qv, n_clu);
                        if (max_w < out_w) {
                            max_w = out_w;
                            max_id = out_id;
                        } else if (max_w == out_w && score[max_id] <= score[out_id]) {
                            max_id = out_id;
                        }
                    }
                    score[cur_id] = max_w + score[max_id];
                    max_out_id[cons_i][cur_id] = max_id;
                }
            }
            for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
                in_id = abg->node[cur_id].in_id[i];
                if (--_out_degree[in_id] == 0) 
                    kdq_push_int(q, in_id);
            }
        }
    }
    abpoa_set_hb_cons(abg, max_out_id, n_clu, clu_read_ids, src_id, sink_id, abc);

    free(score); free(_out_degree);
    for (i = 0; i < n_clu; ++i) free(max_out_id[i]); free(max_out_id);
}

void abpoa_most_freqent(abpoa_graph_t *abg, abpoa_para_t *abpt, int src_id, int sink_id, int *out_degree, int n_clu, int read_ids_n, uint64_t **clu_read_ids, abpoa_cons_t *abc) {
    int use_span_read_count = abpt->sub_aln;
    abpoa_set_msa_rank(abg, src_id, sink_id);
    int i, cons_i, msa_l = abg->node_id_to_msa_rank[sink_id]-1;
    int ***row_column_weight = (int***)_err_malloc(n_clu * sizeof(int**));
    int **msa_node_id = (int**)_err_malloc(msa_l * sizeof(int*));
    // init row_column_weight
    for (cons_i = 0; cons_i < n_clu; ++cons_i) {
        row_column_weight[cons_i] = (int**)_err_malloc(msa_l * sizeof(int*));
        for (i = 0; i < msa_l; ++i) {
            row_column_weight[cons_i][i] = (int*)_err_calloc(abpt->m, sizeof(int));
            row_column_weight[cons_i][i][abpt->m-1] = abc->clu_n_seq[cons_i];
        }
    }
    for (i = 0; i < msa_l; ++i) msa_node_id[i] = (int*)_err_calloc(abpt->m, sizeof(int));
    abc->n_cons = n_clu;
    if (n_clu == 1) abc->clu_n_seq[0] = abc->n_seq;
    else {
        for (cons_i = 0; cons_i < n_clu; ++cons_i) {
            abc->clu_n_seq[cons_i] = get_read_cnt(clu_read_ids[cons_i], read_ids_n);
            set_clu_read_ids(abc, clu_read_ids, cons_i, abc->n_seq);
        }
    }
    // no quality weight used for now, only read count
    abpoa_set_row_column_weight(abg, n_clu, abpt->m, row_column_weight, clu_read_ids, msa_node_id);
    abpoa_set_major_voting_cons(abg, use_span_read_count, abpt->m, row_column_weight, msa_node_id, src_id, sink_id, msa_l, abc);
    for (cons_i = 0; cons_i < n_clu; ++cons_i) {
        for (i = 0; i < msa_l; ++i) {
            free(row_column_weight[cons_i][i]);
        } free(row_column_weight[cons_i]);
    }
    for (i = 0; i < msa_l; ++i) {
        free(msa_node_id[i]);
    } 
    free(row_column_weight); free(msa_node_id);
}

void abpoa_output_fx_consensus(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp) {
    if (out_fp == NULL) return;
    int cons_i, j;
    abpoa_cons_t *abc = ab->abc;
    for (cons_i = 0; cons_i < abc->n_cons; ++cons_i) {
        if (abpt->out_fq) fprintf(out_fp, "@Consensus_sequence");
        else fprintf(out_fp, ">Consensus_sequence");
        if (abc->n_cons > 1) {
            fprintf(out_fp, "_%d ", cons_i+1); // cons_id
            for (j = 0; j < abc->clu_n_seq[cons_i]; ++j) { // cluter read_id
                if (j != 0) fprintf(out_fp, ",");
                fprintf(out_fp, "%d", abc->clu_read_ids[cons_i][j]);
            }
        }
        fprintf(out_fp, "\n");
        for (j = 0; j < abc->cons_len[cons_i]; ++j) {
            fprintf(out_fp, "%c", ab_char256_table[abc->cons_base[cons_i][j]]);
        } fprintf(out_fp, "\n");
        if (abpt->out_fq) {
            fprintf(out_fp, "+Consensus_sequence");
            if (abc->n_cons > 1) {
                fprintf(out_fp, "_%d ", cons_i+1); // cons_id
                for (j = 0; j < abc->clu_n_seq[cons_i]; ++j) { // cluter read_id
                    if (j != 0) fprintf(out_fp, ",");
                    fprintf(out_fp, "%d", abc->clu_read_ids[cons_i][j]);
                }
            }
            fprintf(out_fp, "\n");
            for (j = 0; j < abc->cons_len[cons_i]; ++j) {
                fprintf(out_fp, "%c", abc->cons_phred_score[cons_i][j]);
            } fprintf(out_fp, "\n");
        }
    }
}

abpoa_cons_t *abpoa_allocate_cons(abpoa_cons_t *abc, int n_node, int n_seq, int n_cons) {
    int i;
    abc->n_cons = n_cons, abc->n_seq = n_seq;
    abc->clu_n_seq = (int*)_err_calloc(n_cons, sizeof(int));
    abc->cons_len = (int*)_err_calloc(n_cons, sizeof(int));
    abc->cons_node_ids = (int**)_err_malloc(n_cons * sizeof(int*));
    abc->cons_base = (uint8_t**)_err_malloc(n_cons * sizeof(uint8_t*));
    abc->cons_cov = (int**)_err_malloc(n_cons * sizeof(int*));
    abc->clu_read_ids = (int**)_err_malloc(n_cons * sizeof(int*));
    abc->cons_phred_score = (int**)_err_malloc(n_cons * sizeof(int*));
    for (i = 0; i < n_cons; ++i) {
        abc->cons_node_ids[i] = (int*)_err_malloc(n_node * sizeof(int));
        abc->cons_base[i] = (uint8_t*)_err_malloc(n_node * sizeof(uint8_t));
        abc->cons_cov[i] = (int*)_err_malloc(n_node * sizeof(int));
        abc->clu_read_ids[i] = (int*)_err_malloc(n_seq * sizeof(int));
        abc->cons_phred_score[i] = (int*)_err_malloc(n_node * sizeof(int));
    }
    return abc;
}

int abpoa_collect_cand_het_profile(uint8_t **msa, int msa_l, int n_seq, int m, int min_het, cand_het_t *cand_hets, read_het_profile_t *p, int verbose) {
    int n_het_pos = 0;
    int i, j, k;
    int min_hom = n_seq - min_het; // het >= min_het && <= min_hom
    int *depth = (int*)_err_malloc((m+1) * sizeof(int));

    for (i = 0; i < n_seq; ++i) {
        p[i].read_id = i;
        p[i].start_het_idx = -1; p[i].end_het_idx = -2;
        p[i].alleles = (int*)malloc(msa_l * sizeof(int));
        for (j = 0; j < msa_l; ++j) p[i].alleles[j] = -1; // init as unused
    }

    // print msa
    // for (int i = 0; i < n_seq; ++i) {
    //     fprintf(stderr, ">%d\n", i);
    //     for (int j = 0; j < msa_l; ++j) {
    //         fprintf(stderr, "%c", "ACGTN-"[msa[i][j]]);
    //     } fprintf(stderr, "\n");
    // }

    for (i = 0; i < msa_l; ++i) {
        // if (i == 252)
            // fprintf(stderr, "debug\n");
        memset(depth, 0, (m+1) * sizeof(int));
        for (j = 0; j < n_seq; ++j) {
            depth[msa[j][i]]++;
        }
        int max_base = -1, max_c = 0, sec_base = -1, sec_c = 0;
        for (j = 0; j < m+1; ++j) {
            if (depth[j] > max_c) {
                sec_base = max_base;
                sec_c = max_c;
                max_c = depth[j];
                max_base = j;
            } else if (depth[j] > sec_c) {
                sec_c = depth[j];
                sec_base = j;
            }
        }
        if (max_c >= min_het && max_c <= min_hom && sec_c >= min_het && sec_c <= min_hom) { // X
            cand_hets[n_het_pos].pos = i;
            cand_hets[n_het_pos].phase_set = -1;
            if (max_base == m || sec_base == m) cand_hets[n_het_pos].var_type = 1; // ID
            else cand_hets[n_het_pos].var_type = 0; // X
            cand_hets[n_het_pos].n_depth = n_seq; // XXX
            cand_hets[n_het_pos].n_uniq_alles = 2;
            cand_hets[n_het_pos].alle_covs = (int*)malloc(2 * sizeof(int));
            cand_hets[n_het_pos].alle_covs[0] = max_c; cand_hets[n_het_pos].alle_covs[1] = sec_c;
            cand_hets[n_het_pos].alle_bases = (uint8_t*)malloc(2 * sizeof(uint8_t));
            cand_hets[n_het_pos].alle_bases[0] = max_base; cand_hets[n_het_pos].alle_bases[1] = sec_base;
            n_het_pos++;
        } else continue;
    }
    free(depth);
    for (i = 0; i < n_het_pos; ++i) {
        int pos = cand_hets[i].pos;
        if (verbose >= 2) fprintf(stderr, "het pos: %d, %d (%d) %d (%d)\n", pos, cand_hets[i].alle_bases[0], cand_hets[i].alle_covs[0], cand_hets[i].alle_bases[1], cand_hets[i].alle_covs[1]);
        cand_hets[i].hap_to_alle_profile = NULL;
        cand_hets[i].alle_to_hap = NULL;
        cand_hets[i].hap_to_cons_alle = NULL;
        for (j = 0; j < n_seq; ++j) { // XXX partial aligned reads
            int allele_i = -1;
            for (k = 0; k < cand_hets[i].n_uniq_alles; ++k) {
                if (msa[j][pos] == cand_hets[i].alle_bases[k]) {
                    allele_i = k;
                    break;
                }
            }
            if (allele_i == -1) continue;
            if (p[j].start_het_idx == -1) p[j].start_het_idx = i;
            p[j].end_het_idx = i;
            p[j].alleles[i-p[j].start_het_idx] = allele_i;
        }
    }
    return n_het_pos;
}

int abpoa_collect_cand_het_pos(uint8_t **msa, int msa_l, int n_seq, int m, int min_het, cand_het_pos_t *cand_het_pos, int verbose) {
    int n_het_pos = 0;
    int i, j;
    int min_hom = n_seq - min_het; // het >= min_het && <= min_hom
    int *depth = (int*)_err_malloc((m+1) * sizeof(int));

    // print msa
    // for (int i = 0; i < n_seq; ++i) {
    //     fprintf(stderr, ">%d\n", i);
    //     for (int j = 0; j < msa_l; ++j) {
    //         fprintf(stderr, "%c", "ACGTN-"[msa[i][j]]);
    //     } fprintf(stderr, "\n");
    // }

    for (i = 0; i < msa_l; ++i) {
        // if (i == 252)
            // fprintf(stderr, "debug\n");
        memset(depth, 0, (m+1) * sizeof(int));
        for (j = 0; j < n_seq; ++j) {
            depth[msa[j][i]]++;
        }
        int max_base = -1, max_c = 0, sec_base = -1, sec_c = 0;
        for (j = 0; j < m+1; ++j) {
            if (depth[j] > max_c) {
                sec_base = max_base;
                sec_c = max_c;
                max_c = depth[j];
                max_base = j;
            } else if (depth[j] > sec_c) {
                sec_c = depth[j];
                sec_base = j;
            }
        }
        if (max_c >= min_het && max_c <= min_hom && sec_c >= min_het && sec_c <= min_hom) { // X
            if (verbose >= 2)
                fprintf(stderr, "het pos: %d, %d (%d) %d (%d)\n", i, max_base, max_c, sec_base, sec_c);
            cand_het_pos[n_het_pos].pos = i;
            if (max_base == m || sec_base == m) cand_het_pos[n_het_pos].var_type = 1; // ID
            else cand_het_pos[n_het_pos].var_type = 0; // X
            cand_het_pos[n_het_pos].depth = n_seq; // XXX
            cand_het_pos[n_het_pos].n_uniq_alles = 2;
            cand_het_pos[n_het_pos].alle_bases = (uint8_t*)malloc(2 * sizeof(uint8_t));
            cand_het_pos[n_het_pos].alle_bases[0] = max_base; cand_het_pos[n_het_pos].alle_bases[1] = sec_base;
            n_het_pos++;
        } else continue;
    }
    free(depth);
    return n_het_pos;
}

int abpoa_collect_max_cov_allele(cand_het_t *het) {
    int max_cov = 0, max_cov_alle_i = -1;
    for (int i = 0; i < het->n_uniq_alles; ++i) {
        if (het->alle_covs[i] > max_cov) {
            max_cov = het->alle_covs[i]; max_cov_alle_i = i;
        }
    }
    return max_cov_alle_i;
}

void abpoa_het_init_hap_profile(cand_het_t *hets, int n_cand_hets) {
    for (int het_i = 0; het_i < n_cand_hets; ++het_i) {
        cand_het_t *het = hets+het_i;
        if (het->hap_to_alle_profile == NULL) {
            het->alle_to_hap = (uint8_t*)calloc(het->n_uniq_alles, sizeof(uint8_t)); // +1: minor_alt_allele
            het->hap_to_alle_profile = (int**)malloc(3 * sizeof(int*));
            for (int i = 0; i <= 2; ++i) het->hap_to_alle_profile[i] = (int*)calloc(het->n_uniq_alles, sizeof(int));
            het->hap_to_cons_alle = (int*)malloc(3 * sizeof(int));
            het->hap_to_cons_alle[0] = abpoa_collect_max_cov_allele(het);
            for (int j = 1; j <= 2; ++j) {
                het->hap_to_cons_alle[j] = -1;
            }
        } else {
            memset(het->hap_to_alle_profile[0], 0, het->n_uniq_alles * sizeof(int));
            for (int j = 1; j <= 2; ++j) {
                memset(het->hap_to_alle_profile[j], 0, het->n_uniq_alles * sizeof(int));
            }
            het->hap_to_cons_alle[0] = abpoa_collect_max_cov_allele(het);
            for (int j = 1; j <= 2; ++j) {
                het->hap_to_cons_alle[j] = -1;
            }
        }

        if (het->hap_to_alle_profile == NULL) {
            het->alle_to_hap = (uint8_t*)calloc(1, sizeof(uint8_t)); // +1: minor_alt_allele
            het->hap_to_alle_profile = (int**)malloc(3 * sizeof(int*));
            for (int i = 0; i <= 2; ++i) het->hap_to_alle_profile[i] = (int*)calloc(1, sizeof(int));
            het->hap_to_cons_alle = (int*)malloc(3 * sizeof(int));
            het->hap_to_cons_alle[0] = 0;
            for (int j = 1; j <= 2; ++j) {
                het->hap_to_cons_alle[j] = -1;
            }
        } else {
            memset(het->hap_to_alle_profile[0], 0, 1* sizeof(int));
            for (int j = 1; j <= 2; ++j) {
                memset(het->hap_to_alle_profile[j], 0, 1 * sizeof(int));
            }
            het->hap_to_cons_alle[0] = 0;
            for (int j = 1; j <= 2; ++j) {
                het->hap_to_cons_alle[j] = -1;
            }
        }
    }
}

int *sort_cand_hets(cand_het_t *cand_hets, int n_cand_hets) {
    int *sorted_cand_het_i = (int*)malloc(n_cand_hets * sizeof(int));
    for (int i = 0; i < n_cand_hets; ++i) {
        sorted_cand_het_i[i] = i;
    }
    // sort var by type: clean SNP -> clean indel -> noisy SNP -> noisy INDEL
    for (int i = 0; i < n_cand_hets-1; ++i) {
        for (int j = i+1; j < n_cand_hets; ++j) {
            if (cand_hets[i].var_type > cand_hets[j].var_type) {
                int tmp = sorted_cand_het_i[i]; sorted_cand_het_i[i] = sorted_cand_het_i[j]; sorted_cand_het_i[j] = tmp;
            }
        }
    }
    return sorted_cand_het_i;
}

int abpoa_assign_het_init_hap(cand_het_t *het, int verbose) {
    if (verbose >= 2) fprintf(stderr, "Init Het-Var hap: %d %c\n", het->pos, "XG"[het->var_type]);
    if (het->var_type == 0) het->phase_set = het->pos; // potential start of a PhaseSet
    else het->phase_set = het->pos - 1;
    int hap1_alle_i = 0, hap2_alle_i = 1;
    het->alle_to_hap[hap1_alle_i] = 1; het->alle_to_hap[hap2_alle_i] = 2;
    return het->phase_set;
}

// assign haplotype to a SNP based on read SNP profiles
// 1. pick the most common haplotype and corresponding most common base 
// 2. assign most common base to the most common haplotype, and other bases to the other haplotype
int abpoa_assign_het_hap_based_on_pre_reads1(cand_het_t *het, int verbose) {
    int first_hap=0, first_hap_cnt=0, first_hap_alle_i=-1;
    int sec_hap=0, sec_hap_cnt=0, sec_hap_alle_i=-1;

    for (int j = 1; j <= 2; ++j) {
        for (int i = 0; i < het->n_uniq_alles; ++i) { // ref + alt_allele; minor_alt_allele is not considered
            if (het->hap_to_alle_profile[j][i] > first_hap_cnt) {
                sec_hap_cnt = first_hap_cnt; sec_hap = first_hap; sec_hap_alle_i = first_hap_alle_i;
                first_hap_cnt = het->hap_to_alle_profile[j][i]; first_hap = j; first_hap_alle_i = i;
            } else if (het->hap_to_alle_profile[j][i] > sec_hap_cnt) {
                sec_hap_cnt = het->hap_to_alle_profile[j][i]; sec_hap = j; sec_hap_alle_i = i;
            }
        }
    }
    if (sec_hap == 0) {
        if (first_hap == 1) sec_hap = 2; else sec_hap = 1;
        // set the most common bases other than first_hap_base to sec_hap
        int sec_hap_alle_cov = 0;
        for (int i = 0; i < het->n_uniq_alles; ++i) {
            if (het->alle_covs[i] > sec_hap_alle_cov && i != first_hap_alle_i) {
                sec_hap_alle_i = i; sec_hap_alle_cov = het->alle_covs[i];
            }
        }
    }
    if (first_hap == sec_hap) {
        if (verbose >= 2) fprintf(stderr, "Var: %d, %c, first_hap: %d (%d: %d), sec_hap: %d (%d: %d)\n", het->pos, "XG"[het->var_type], first_hap, first_hap_alle_i, first_hap_cnt, sec_hap, sec_hap_alle_i, sec_hap_cnt);
        abpoa_assign_het_init_hap(het, verbose);
    // } else if (first_hap_alle_i == sec_hap_alle_i) { // homozygous
    } else {
        for (int i = 0; i < het->n_uniq_alles; ++i) {
            if (i == first_hap_alle_i) het->alle_to_hap[i] = first_hap;
            else if (i == sec_hap_alle_i) het->alle_to_hap[i] = sec_hap;
            else het->alle_to_hap[i] = 0;
        }
    }
    return het->phase_set;
} 

void abpoa_het_init_hap_cons_alle0(cand_het_t *het, int min_alt_dp) {
    // select the most common allele as the consensus allele, based on hap_to_alle_profile
    for (int hap = 1; hap <= 2; ++hap) {
        int max_cov = 0, max_cov_alle_i = -1;
        for (int i = 0; i < het->n_uniq_alles; ++i) {
            if (het->hap_to_alle_profile[hap][i] > max_cov) {
                max_cov = het->hap_to_alle_profile[hap][i];
                if (max_cov >= min_alt_dp) max_cov_alle_i = i;
            }
        }
        het->hap_to_cons_alle[hap] = max_cov_alle_i;
    }
}

int abpoa_het_hap_profile_cov(cand_het_t *het) {
    int cov = 0;
    for (int j = 1; j <= 2; ++j) {
        for (int i = 0; i < het->n_uniq_alles; ++i) {
            cov += het->hap_to_alle_profile[j][i];
        }
    }
    return cov;
}

int abpoa_assign_het_hap_based_on_pre_reads(cand_het_t *het, int min_dp, int verbose) {
    if (abpoa_het_hap_profile_cov(het) < min_dp) 
        return abpoa_assign_het_init_hap(het, verbose);
    else {
        return abpoa_assign_het_hap_based_on_pre_reads1(het, verbose);
    }
}


// after a read is assigned with hap, update hap of all other SNPs covered by this read
// including homozygous and heterozygous SNPs
void abpoa_update_het_hap_profile_based_on_aln_hap(int hap, int phase_set, cand_het_t *hets, read_het_profile_t *p, int read_i) {
    int start_het_idx = p[read_i].start_het_idx, end_het_idx = p[read_i].end_het_idx;
    for (int het_i = start_het_idx; het_i <= end_het_idx; ++het_i) {
        int read_het_idx = het_i - start_het_idx;
        // if (p[read_i].var_is_used[read_var_idx] == 0) continue;
        int allele_i = p[read_i].alleles[read_het_idx];
        if (allele_i == -1) continue;
        hets[het_i].hap_to_alle_profile[hap][allele_i] += 1;
        // update HOM var's phase_set as well ??
        if (hets[het_i].phase_set == -1 || phase_set <= hets[het_i].pos)
            hets[het_i].phase_set = phase_set;
    }
}

int abpoa_collect_tmp_hap_cons_allele_by_deduct_read(cand_het_t *het, int hap, int allele_i, int *tmp_hap_to_cons_alle, int verbose) {
    for (int i = 1; i <= 2; ++i) {
        tmp_hap_to_cons_alle[i] = -1;
        if (i != hap || het->hap_to_cons_alle[i] != allele_i) { // no change
            tmp_hap_to_cons_alle[i] = het->hap_to_cons_alle[i];
        } else { // i==hap && var->hap_to_cons_alle[i] == allele_i
            int cov, max_cov = 0, max_cov_alle_i = -1;
            for (int j = 0; j < het->n_uniq_alles; ++j) {
                cov = het->hap_to_alle_profile[i][j];
                if (j == allele_i) cov -= 1;
                if (cov > max_cov) {
                    max_cov = cov; max_cov_alle_i = j;
                }
            }
            if (max_cov_alle_i == -1) {
                if (verbose >= 2) fprintf(stderr, "No HAP %d allele in SNP: %d, %c\n", i, het->pos, "XG"[het->var_type]);
            }
            tmp_hap_to_cons_alle[i] = max_cov_alle_i;
        }
    }
    return 0;
}

// update hap_cons_profile based on hap_to_base_profile
// update haplotype for a read based on SNP profiles of all other overlapping reads
// input: target read, SNP profiles of all reads
// output: updated haplotype of the target read
int abpoa_update_het_aln_hap1(int target_read_i, int cur_hap, read_het_profile_t *p, cand_het_t *cand_hets, int verbose) {
    int start_het_idx = p[target_read_i].start_het_idx, end_het_idx = p[target_read_i].end_het_idx;
    // deduct target read from hap_to_alle_profile, then compare target read's var profile with hap_cons_alle
    int *hap_match_cnt = (int*)calloc(3, sizeof(int));
    int *tmp_hap_to_cons_alle = (int*)malloc(3 * sizeof(int));

    for (int het_i = start_het_idx; het_i <= end_het_idx; ++het_i) {
        int read_het_idx = het_i - start_het_idx;
        cand_het_t *het = cand_hets+het_i;
        int var_weight = 1; // XXX
        if (het->var_type == 0) var_weight = 2;

        int allele_i = p[target_read_i].alleles[read_het_idx];
        if (allele_i == -1) continue;
        abpoa_collect_tmp_hap_cons_allele_by_deduct_read(het, cur_hap, allele_i, tmp_hap_to_cons_alle, verbose);
        for (int i = 1; i <= 2; ++i) {
            if (tmp_hap_to_cons_alle[i] == allele_i) {
                hap_match_cnt[i] += var_weight;
            }
        }
    }
    int max_cnt=0, max_hap=0, sec_cnt=0;
    for (int i = 1; i <= 2; ++i) {
        if (hap_match_cnt[i] > max_cnt) {
            sec_cnt = max_cnt;
            max_cnt = hap_match_cnt[i]; max_hap = i;
        } else if (hap_match_cnt[i] > sec_cnt) {
            sec_cnt = hap_match_cnt[i];
        }
    }
    free(hap_match_cnt); free(tmp_hap_to_cons_alle);
    if (max_cnt == 0) {
        fprintf(stderr, "Read %d max_cnt == 0\n", target_read_i);
        return 0; // unknown
    } else if (max_cnt == sec_cnt) {
        fprintf(stderr, "Read %d max_cnt == sec_cnt\n", target_read_i);
        max_hap = cur_hap;
    }
    return max_hap;
}

int abpoa_update_het_hap_profile_based_on_changed_hap(int new_hap, int old_hap, cand_het_t *cand_hets, int min_alt_dp, read_het_profile_t *p, int read_i, int verbose) {
    int start_het_idx = p[read_i].start_het_idx, end_het_idx = p[read_i].end_het_idx;
    for (int het_i = start_het_idx; het_i <= end_het_idx; ++het_i) {
        int read_het_idx = het_i - start_het_idx;
        // if (p[read_i].var_is_used[read_var_idx] == 0) continue;
        int allele_i = p[read_i].alleles[read_het_idx];
        if (allele_i == -1) continue;
        cand_het_t *het = cand_hets+het_i;
        // only update Het-Var
        if (verbose >= 2) fprintf(stderr, "pos: %d, old_hap: %d, new_hap: %d, var: %d\n", het->pos, old_hap, new_hap, allele_i);
        het->hap_to_alle_profile[old_hap][allele_i] -= 1;
        het->hap_to_alle_profile[new_hap][allele_i] += 1;
        abpoa_het_init_hap_cons_alle0(het, min_alt_dp);
    }
    return 0;
}

// TODO: > 2 clusters
void abpoa_clu_reads_based_on_het_pos(int n_het_pos, cand_het_t *cand_hets, read_het_profile_t *p, uint64_t ***clu_read_ids, int n_seq, int verbose) {
    int *haps = (int*)_err_calloc(n_seq, sizeof(int));
    (*clu_read_ids) = (uint64_t**)_err_malloc(sizeof(uint64_t*) * 2);
    int read_id_n = (n_seq-1)/64+1;
    for (int i = 0; i < 2; ++i) {
        (*clu_read_ids)[i] = (uint64_t*)_err_calloc(read_id_n, sizeof(uint64_t));
    }
    abpoa_het_init_hap_profile(cand_hets, n_het_pos);
    int *sorted_cand_het_i = sort_cand_hets(cand_hets, n_het_pos);
    for (int _het_i = 0; _het_i < n_het_pos; ++_het_i) {
        int het_i = sorted_cand_het_i[_het_i];
        cand_het_t *het = cand_hets+het_i;
        // here we highly rely on reads that were assigned with HAPs previously, >= 2 reads is enough
        int phase_set = abpoa_assign_het_hap_based_on_pre_reads(het, 2, verbose); //opt->min_dp); // update alle_to_hap
        for (int read_i = 0; read_i < n_seq; ++read_i) {
            if (p[read_i].start_het_idx == -1) continue;
            int read_het_idx = het_i - p[read_i].start_het_idx;
            if (haps[read_i] == 0) {
                int het_alle_i = p[read_i].alleles[read_het_idx];
                if (het_alle_i == -1) continue;
                int hap = het->alle_to_hap[het_alle_i];
                if (hap != 0) {
                    // first time assign hap to the read (update bam_haps)
                    haps[read_i] = hap;
                    if (verbose >= 2) fprintf(stderr, "read: %d, cur_var: %d, %c, alle: %d, hap: %d\n", read_i, het->pos, "XG"[het->var_type], het_alle_i, hap);
                    // update hap_to_alle_profile for all Vars covered by this read, based on its assigned haplotype
                    // udpated profile will then be used for following Vars (assign_var_hap_based_on_pre_reads)
                    abpoa_update_het_hap_profile_based_on_aln_hap(hap, phase_set, cand_hets, p, read_i);
                }
            }
        }
        abpoa_het_init_hap_cons_alle0(het, 2); // update hap_to_cons_alle
    } // after first round, 
    // bam_haps/hap_to_alle_profile/hap_to_cons_alle are upToDate and will be used in the following rounds
    // alle_to_hap will not be used (may be NOT upToDate)
    // 2nd loop: read-wise iterative loop
    int changed_hap, max_iter = 10, i_iter=0;
    while (i_iter++ < max_iter) {
        if (verbose >= 2) fprintf(stderr, "iter: %d\n", i_iter);
        changed_hap = 0;
        // read-wise loop
        // re-calculate read-wise haplotype (Var-wise 1. Hap, 2. Allele, 3. Read, 4. AlleleCons are all up-to-date)
        for (int read_i = 0; read_i < n_seq; ++read_i) {
            if (haps[read_i] == 0) continue;
            int cur_hap = haps[read_i];
            // XXX TODO: potential local optima
            int new_hap = abpoa_update_het_aln_hap1(read_i, cur_hap, p, cand_hets, verbose);
            if (new_hap != cur_hap) { // update bam_haps, hap_to_base_profile, hap_to_cons_base
                if (verbose >= 2) {
                    fprintf(stderr, "read %d:\t", read_i);
                    fprintf(stderr, "\t\t cur_hap: %d, new_hap: %d\n", cur_hap, new_hap);
                }
                changed_hap = 1;
                haps[read_i] = new_hap; // update intermediately
                // bam_aux_append(chunk->reads[read_i], "XT", 'i', 4, (uint8_t*)&(chunk->haps[read_i]));
                abpoa_update_het_hap_profile_based_on_changed_hap(new_hap, cur_hap, cand_hets, 2, p, read_i, verbose);
            }
        } if (changed_hap == 0) break;
    }
    // assign read to cluster based on haplotype
    for (int read_i = 0; read_i < n_seq; ++read_i) {
        if (haps[read_i] == 0) continue;
        int hap = haps[read_i];
        (*clu_read_ids)[hap-1][read_i/64] |= (1ULL << (read_i & 0x3f));
    }
    free(sorted_cand_het_i); free(haps);
}

int abpoa_collect_msa_dis1(uint8_t **msa, int msa_l, int n_seq, int m, cand_het_pos_t *cand_het_pos, int n_het_pos, int i, int j) {
    int dis = 0;
    for (int k = 0; k < n_het_pos; ++k) {
        int pos = cand_het_pos[k].pos;
        int var_weight = 1; // Gap
        if (cand_het_pos[k].var_type == 0) var_weight = 2; // SNP
        uint8_t base1 = msa[i][pos], base2 = msa[j][pos];
        uint8_t alle_base1 = cand_het_pos[k].alle_bases[0], alle_base2 = cand_het_pos[k].alle_bases[1];
        if (base1 != alle_base1 && base1 != alle_base2) continue;
        if (base2 != alle_base1 && base2 != alle_base2) continue;
        if (base1 != base2) dis += var_weight;
    }
    return dis;
}

int **abpoa_collect_msa_dis_matrix(uint8_t **msa, int msa_l, int n_seq, int m, cand_het_pos_t *cand_het_pos, int n_het_pos, int verbose) {
    int i, j, **dis_matrix = (int**)malloc(n_seq * sizeof(int*));
    for (i = 0; i < n_seq; ++i) {
        dis_matrix[i] = (int*)calloc(n_seq, sizeof(int));
    }

    for (i = 0; i < n_seq; ++i) {
        for (j = i+1; j < n_seq; ++j) {
            dis_matrix[i][j] = abpoa_collect_msa_dis1(msa, msa_l, n_seq, m, cand_het_pos, n_het_pos, i, j);
            dis_matrix[j][i] = dis_matrix[i][j];
        }
    }
    return dis_matrix;
}

// XXX TODO: > 2 clusters
int abpoa_init_kmedoids(int **dis_matrix, int n_seq, int max_n_cons, int **medoids) {
    int max_dis = 0, med1 = -1, med2 = -1;
    for (int i = 0; i < n_seq-1; ++i) {
        for (int j = i+1; j < n_seq; ++j) {
            if (dis_matrix[i][j] > max_dis) {
                max_dis = dis_matrix[i][j];
                med1 = i; med2 = j;
            }
        }
    }
    if (max_dis > 0) {
        (*medoids) = (int*)malloc(2 * sizeof(int));
        (*medoids)[0] = med1; (*medoids)[1] = med2;
        return 0;
    } else return -1;
}

// within each cluster, find the medoid that has the smallest sum of distances to all other members
void abpoa_collect_kmedoids0(int **dis_matrix, int n_seq, int max_n_cons, int **clu_reads, int *n_clu_seqs, int *medoids) {
    for (int i = 0; i < max_n_cons; ++i) {
        int min_sum_dis = INT_MAX, min_sum_dis_read_i = -1;
        for (int j = 0; j < n_clu_seqs[i]; ++j) {
            int sum_dis = 0, read_i = clu_reads[i][j];
            for (int k = 0; k < n_clu_seqs[i]; ++k) {
                if (j == k) continue;
                int read_j = clu_reads[i][k];
                sum_dis += dis_matrix[read_i][read_j];
            }
            if (sum_dis < min_sum_dis) {
                min_sum_dis = sum_dis;
                min_sum_dis_read_i = read_i;
            }
        }
        if (min_sum_dis_read_i != -1) medoids[i] = min_sum_dis_read_i;
    }
    // sort medoids, from small to large
    for (int i = 0; i < max_n_cons-1; ++i) {
        for (int j = i+1; j < max_n_cons; ++j) {
            if (medoids[i] > medoids[j]) {
                int tmp = medoids[i]; medoids[i] = medoids[j]; medoids[j] = tmp;
            }
        }
    }
}

int abpoa_update_kmedoids(int **dis_matrix, int n_seq, int max_n_cons, int **medoids, int **clu_reads, int *n_clu_seqs, int verbose) {
    int *new_medoids = (int*)malloc(max_n_cons * sizeof(int));
    for (int i = 0; i < max_n_cons; ++i) new_medoids[i] = -1;
    n_clu_seqs[0] = 0; n_clu_seqs[1] = 0;
    for (int i = 0; i < n_seq; ++i) {
        int min_dis = INT_MAX, min_clu = -1, tied = 0;
        for (int j = 0; j < max_n_cons; ++j) {
            if (dis_matrix[i][(*medoids)[j]] < min_dis) {
                min_dis = dis_matrix[i][(*medoids)[j]];
                min_clu = j;
                tied = 0;
            } else if (dis_matrix[i][(*medoids)[j]] == min_dis) {
                tied = 1;
            }
        }
        if (tied != 0 || min_clu == -1) continue;
        clu_reads[min_clu][n_clu_seqs[min_clu]++] = i;
    }
    abpoa_collect_kmedoids0(dis_matrix, n_seq, max_n_cons, clu_reads, n_clu_seqs, new_medoids);
    int changed = 0;
    for (int i = 0; i < max_n_cons; ++i) {
        if (new_medoids[i] == -1) { // no reads in this cluster
            changed = 0; break;
        }
        if (new_medoids[i] != (*medoids)[i]) changed = 1;
    } 
    free(*medoids);
    *medoids = new_medoids;
    return changed;
}

int abpoa_clu_reads_kmedoids(int **dis_matrix, int n_seq, int min_het, int max_n_cons, uint64_t ***clu_read_ids, int verbose) {
    // init
    int *medoids = NULL;
    if (abpoa_init_kmedoids(dis_matrix, n_seq, max_n_cons, &medoids) == -1) return 0;

    // iterative update clusters until convergence
    int **clu_reads = (int**)_err_malloc(sizeof(int*) * max_n_cons);
    int *n_clu_seqs = (int*)_err_malloc(sizeof(int) * max_n_cons);
    for (int i = 0; i < max_n_cons; ++i) {
        clu_reads[i] = (int*)_err_calloc(n_seq, sizeof(int));
    }
    int iter = 0;
    while (1) {
        int changed = abpoa_update_kmedoids(dis_matrix, n_seq, max_n_cons, &medoids, clu_reads, n_clu_seqs, verbose);
        iter++;
        if (changed == 0 || iter >= 10) break;
    }
    int n_clu = 0, n_clustered_reads = 0;
    for (int i = 0; i < max_n_cons; ++i) {
        if (n_clu_seqs[i] >= min_het) n_clu++;
        n_clustered_reads += n_clu_seqs[i];
    }
    if (n_clu != max_n_cons || n_clustered_reads < ceil(n_seq * 0.8)) n_clu = 1;
    else {
        (*clu_read_ids) = (uint64_t**)_err_malloc(sizeof(uint64_t*) * n_clu);
        int read_id_n = (n_seq-1)/64+1;
        for (int i = 0; i < n_clu; ++i) {
            (*clu_read_ids)[i] = (uint64_t*)_err_calloc(read_id_n, sizeof(uint64_t));
            for (int j = 0; j < n_clu_seqs[i]; ++j) {
                int read_i = clu_reads[i][j];
                (*clu_read_ids)[i][read_i/64] |= (1ULL << (read_i & 0x3f));
            }
        }
    }
    for (int i = 0; i < max_n_cons; ++i) free(clu_reads[i]);
    free(clu_reads); free(n_clu_seqs); free(medoids);
    return n_clu;
}

int abpoa_multip_read_clu_kmedoids(abpoa_graph_t *abg, abpoa_para_t *abpt, int src_id, int sink_id, int n_seq, int m, int max_n_cons, double min_freq, uint64_t ***clu_read_ids, int verbose) {
    abpoa_set_msa_rank(abg, src_id, sink_id);
    int i, n_clu;
    uint8_t **msa = (uint8_t**)_err_malloc(sizeof(uint8_t*) * n_seq);
    int msa_l = abpoa_collect_msa(abg, abpt, msa, n_seq); // abg->node_id_to_msa_rank[sink_id]-1, 
    int min_w = MAX_OF_TWO(2, ceil(n_seq * min_freq));
    cand_het_pos_t *cand_het_pos = (cand_het_pos_t*)_err_malloc(msa_l * sizeof(cand_het_pos_t));
    int n_het_pos = abpoa_collect_cand_het_pos(msa, msa_l, n_seq, m, min_w, cand_het_pos, verbose);
    if (n_het_pos < 1) n_clu = 1;
    else {
        int **dis_matrix = abpoa_collect_msa_dis_matrix(msa, msa_l, n_seq, m, cand_het_pos, n_het_pos, verbose);
        n_clu = abpoa_clu_reads_kmedoids(dis_matrix, n_seq, min_w, max_n_cons, clu_read_ids, verbose);
        for (i = 0; i < n_seq; ++i) free(dis_matrix[i]); free(dis_matrix);
    }
    for (int i = 0; i < n_het_pos; ++i) free(cand_het_pos[i].alle_bases); free(cand_het_pos);
    for (int i = 0; i < n_seq; ++i) free(msa[i]); free(msa);
    return n_clu;
}
// XXX: more than 2 clusters
int abpoa_multip_read_clu2(abpoa_graph_t *abg, abpoa_para_t *abpt, int src_id, int sink_id, int n_seq, int m, int max_n_cons, double min_freq, uint64_t ***clu_read_ids, int verbose) {
    abpoa_set_msa_rank(abg, src_id, sink_id);
    int i, j, n_clu;
    uint8_t **msa = (uint8_t**)_err_malloc(sizeof(uint8_t*) * n_seq);
    int msa_l = abpoa_collect_msa(abg, abpt, msa, n_seq); // abg->node_id_to_msa_rank[sink_id]-1, 
    int min_w = MAX_OF_TWO(2, ceil(n_seq * min_freq));
    
    cand_het_t *cand_hets = (cand_het_t*)_err_malloc(msa_l * sizeof(cand_het_t));
    read_het_profile_t *p = (read_het_profile_t*)_err_malloc(n_seq * sizeof(read_het_profile_t));
    int n_het_pos = abpoa_collect_cand_het_profile(msa, msa_l, n_seq, m, min_w, cand_hets, p, verbose);
    for (i = 0; i < n_seq; ++i) free(msa[i]); free(msa);
    
    if (n_het_pos < 1) n_clu = 1;
    else {
        n_clu = 2;
        abpoa_clu_reads_based_on_het_pos(n_het_pos, cand_hets, p, clu_read_ids, n_seq, verbose);
    }
    for (i = 0; i < n_het_pos; ++i) {
        free(cand_hets[i].alle_covs); free(cand_hets[i].alle_bases);
        free(cand_hets[i].alle_to_hap); free(cand_hets[i].hap_to_cons_alle); 
        for (j = 0; j < 3; ++j) free(cand_hets[i].hap_to_alle_profile[j]);
        free(cand_hets[i].hap_to_alle_profile);
    } free(cand_hets);
    for (i = 0; i < n_seq; ++i) {
        free(p[i].alleles);
    } free(p);
    return n_clu;
}

void abpoa_generate_2cons(abpoa_t *ab, abpoa_para_t *abpt) {
    return;
}

// should always do topological sort first, then generate consensus
void abpoa_generate_consensus(abpoa_t *ab, abpoa_para_t *abpt) {
    if (ab->abg->is_called_cons == 1) return;
    abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) return;
    // if (abpt->max_n_cons == 2) {
    //     abpoa_generate_2cons(ab, abpt);
    //     return;
    // }
    int i, *out_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    for (i = 0; i < abg->node_n; ++i) {
        out_degree[i] = abg->node[i].out_edge_n;
    }

    int n_clu, n_seq = ab->abs->n_seq; uint64_t **clu_read_ids=NULL;
    int read_ids_n = (n_seq-1)/64+1;

    n_clu = 1;
    if (abpt->max_n_cons > 1) {
        n_clu = abpoa_multip_read_clu_kmedoids(abg, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, n_seq, abpt->m, abpt->max_n_cons, abpt->min_freq, &clu_read_ids, abpt->verbose);
    }

    abpoa_cons_t *abc = ab->abc;
    abpoa_allocate_cons(abc, abg->node_n, ab->abs->n_seq, n_clu);
    if (abpt->cons_algrm == ABPOA_HB)
        abpoa_heaviest_bundling(abg, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree, n_clu, read_ids_n, clu_read_ids, abc);
    else
        abpoa_most_freqent(abg, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree, n_clu, read_ids_n, clu_read_ids, abc);
    if (n_clu > 1) {
        for (i = 0; i < n_clu; ++i) free(clu_read_ids[i]); free(clu_read_ids);
    }
    abg->is_called_cons = 1; free(out_degree);
}
