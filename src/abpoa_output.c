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

void set_clu_read_ids0(abpoa_cons_t *abc, uint64_t **read_ids, int n_seq) {
    int i;
    for (i = 0; i < n_seq; ++i) {
        abc->clu_read_ids[0][i] = i;
    }
    if (i != abc->clu_n_seq[0])
        err_fatal(__func__, "Error in set cluster read ids. (%d, %d)", i, abc->clu_n_seq[0]);
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
    if (n_clu == 1) {
        abc->clu_n_seq[0] = abc->n_seq;
        set_clu_read_ids0(abc, clu_read_ids, abc->n_seq);
    } else {
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

void abpoa_most_frequent(abpoa_graph_t *abg, abpoa_para_t *abpt, int src_id, int sink_id, int n_clu, int read_ids_n, uint64_t **clu_read_ids, abpoa_cons_t *abc) {
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
    if (n_clu == 1) {
        abc->clu_n_seq[0] = abc->n_seq;
        set_clu_read_ids0(abc, clu_read_ids, abc->n_seq);
    } else {
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
        if (abpt->batch_index > 0) {
            fprintf(out_fp, "_%d", abpt->batch_index); // batch index for file mapping
        }
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
            if (abpt->batch_index > 0) {
                fprintf(out_fp, "_%d", abpt->batch_index); // batch index for file mapping
            }
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

int allele_clu_exist(cand_het_pos_t *cand_het_pos, int n_het_pos, int n_uniq_alles, int *alleles, int *n_clu_reads, int **clu_read_ids) {
    int het_i = -1;
    int i, j, k, exist;
    for (i = n_het_pos-1; i >= 0; --i) {
        if (cand_het_pos[i].n_uniq_alles != n_uniq_alles) continue;
        exist = 1;
        for (j = 0; j < n_uniq_alles; ++j) {
            int allele_x = cand_het_pos[i].alle_bases[j];
            int allele_y = alleles[j];
            if (cand_het_pos[i].n_clu_reads[allele_x] != n_clu_reads[allele_y]) {
                exist = 0; break;
            }
            for (k = 0; k < n_clu_reads[allele_y]; ++k) {
                if (cand_het_pos[i].clu_read_ids[allele_x][k] != clu_read_ids[allele_y][k]) {
                    exist = 0; break;
                }
            }
            if (exist == 0) break;
        }
        if (exist) {
            het_i = i; break;
        }
    }
    return het_i;
}

// allow >2 alleles for each het position
int abpoa_collect_cand_het_pos(uint8_t **msa, int msa_l, int n_seq, int m, int min_het, cand_het_pos_t *cand_het_pos, int *cand_het_priority, int verbose) {
    int n_het_pos = 0;
    int i, j, k;
    // int min_hom = n_seq - min_het; // het >= min_het && <= min_hom
    min_het = MAX_OF_TWO(2, min_het/2);
    int min_hom = n_seq - min_het; // het >= min_het && <= min_hom
    int *depth = (int*)_err_malloc((m+1) * sizeof(int));
    int n_uniq_alles, var_type;
    int *alleles = (int*)malloc((m+1) * sizeof(int)); // 0..m, alleles are sorted by the first appearance in n_seq reads
    int *allele_to_idx = (int*)malloc((m+1) * sizeof(int));
    int *allele_apperance = (int*)malloc((m+1) * sizeof(int)); // 0..m, record the first appearance of each allele in n_seq reads
    int *n_clu_reads = (int*)malloc((m+1) * sizeof(int));
    int **clu_read_ids = (int**)malloc((m+1) * sizeof(int*));
    for (j = 0; j < m+1; ++j) {
        clu_read_ids[j] = (int*)malloc(n_seq * sizeof(int));
    }
    for (i = 0; i < msa_l; ++i) {
        n_uniq_alles = 0; var_type = 0;
        int total_depth = 0;
        memset(depth, 0, (m+1) * sizeof(int));
        memset(n_clu_reads, 0, (m+1) * sizeof(int));
        for (j = 0; j < n_seq; ++j) {
            depth[msa[j][i]]++;
            if (depth[msa[j][i]] == 1) // first appearance
                allele_apperance[msa[j][i]] = j;
        }
        for (j = 0; j < m+1; ++j) {
            if (depth[j] >= min_het && depth[j] <= min_hom) {
                alleles[n_uniq_alles] = j;
                total_depth += depth[j];
                n_uniq_alles++;
                if (j == m) var_type = 1; // indel
            }
        }
        if (n_uniq_alles < 2) continue;
        // sort alleles by the first appearance in n_seq reads
        for (j = 0; j < n_uniq_alles-1; ++j) {
            for (k = j+1; k < n_uniq_alles; ++k) {
                if (allele_apperance[alleles[j]] > allele_apperance[alleles[k]]) {
                    int tmp = alleles[j];
                    alleles[j] = alleles[k];
                    alleles[k] = tmp;
                }
            }
        }
        for (j = 0; j < n_uniq_alles; ++j) allele_to_idx[alleles[j]] = j;

        for (j = 0; j < n_seq; ++j) {
            for (k = 0; k < n_uniq_alles; ++k) {
                if (msa[j][i] == alleles[k]) {
                    clu_read_ids[alleles[k]][n_clu_reads[alleles[k]]++] = j;
                    break;
                }
            }
        }
        int het_i = allele_clu_exist(cand_het_pos, n_het_pos, n_uniq_alles, alleles, n_clu_reads, clu_read_ids);
        if (het_i >= 0) { // exist
            // fprintf(stdout, "het pos: %d, same as %d\n", i, cand_het_pos[het_i].pos);
            cand_het_pos[het_i].count++;
            if (var_type == 0) cand_het_pos[het_i].var_type = 0; // X if any is X
            continue;
        }
        // new het pos
        // fprintf(stdout, "het pos: %d, ", i);
        // for (j = 0; j < n_uniq_alles; ++j) {
        //     fprintf(stdout, "%d (%d): ", alleles[j], n_clu_reads[alleles[j]]);
        //     for (k = 0; k < n_clu_reads[alleles[j]]; ++k) {
        //         if (k != 0) fprintf(stdout, ",");
        //         fprintf(stdout, "%d", clu_read_ids[alleles[j]][k]);
        //     } fprintf(stdout, "; ");
        // } fprintf(stdout, "\n");
        cand_het_pos[n_het_pos].pos = i;
        cand_het_pos[n_het_pos].depth = total_depth; // XXX partial aligned reads???
        cand_het_pos[n_het_pos].var_type = var_type; // ID
        cand_het_pos[n_het_pos].count = 1;
        cand_het_pos[n_het_pos].n_uniq_alles = n_uniq_alles;
        cand_het_pos[n_het_pos].n_clu_reads = (int*)malloc((m+1)* sizeof(int));
        cand_het_pos[n_het_pos].clu_read_ids = (int**)malloc((m+1)* sizeof(int*));
        cand_het_pos[n_het_pos].read_id_to_allele_idx = (int*)malloc(n_seq * sizeof(int));
        for (j = 0; j < n_seq; ++j) cand_het_pos[n_het_pos].read_id_to_allele_idx[j] = -1; // init
        for (j = 0; j < m+1; ++j) {
            cand_het_pos[n_het_pos].n_clu_reads[j] = n_clu_reads[j];
            cand_het_pos[n_het_pos].clu_read_ids[j] = (int*)malloc(n_seq * sizeof(int));
            memcpy(cand_het_pos[n_het_pos].clu_read_ids[j], clu_read_ids[j], n_clu_reads[j] * sizeof(int));
            for (k = 0; k < n_clu_reads[j]; ++k) {
                cand_het_pos[n_het_pos].read_id_to_allele_idx[clu_read_ids[j][k]] = allele_to_idx[j];
            }
        }
        cand_het_pos[n_het_pos].alle_bases = (uint8_t*)malloc(n_uniq_alles * sizeof(uint8_t));
        for (j = 0; j < n_uniq_alles; ++j) {
            cand_het_pos[n_het_pos].alle_bases[j] = alleles[j];
        }
        cand_het_priority[n_het_pos] = n_het_pos;
        n_het_pos++;
    }
    // sort cand_het_pos by count and type: var_type (X > ID), store new priority in cand_het_priority
    for (i = 0; i < n_het_pos; ++i) {
        cand_het_priority[i] = i;
    }
    // sort by count, depth, var_type using bubble sort
    int swapped;
    do {
        swapped = 0;
        for (j = 0; j < n_het_pos-1; ++j) {
            int _j = cand_het_priority[j];
            int _j1 = cand_het_priority[j+1];
            if (cand_het_pos[_j].count < cand_het_pos[_j1].count || (cand_het_pos[_j].count == cand_het_pos[_j1].count && cand_het_pos[_j].depth < cand_het_pos[_j1].depth) ||
                (cand_het_pos[_j].count == cand_het_pos[_j1].count && cand_het_pos[_j].depth == cand_het_pos[_j1].depth && cand_het_pos[_j].var_type > cand_het_pos[_j1].var_type)) {
                // printf("Swap: %d (%d, %d) <-> %d (%d, %d)\n", j, cand_het_pos[_j].count, cand_het_pos[_j].var_type, j+1, cand_het_pos[_j1].count, cand_het_pos[_j1].var_type);
                int tmp = cand_het_priority[j];
                cand_het_priority[j] = cand_het_priority[j+1];
                cand_het_priority[j+1] = tmp;
                swapped = 1;
            }
        }
    } while (swapped);

    // selection sort
    // for (i = 0; i < n_het_pos-1; ++i) {
    //     int _i = cand_het_priority[i];
    //     for (j = i+1; j < n_het_pos; ++j) {
    //         int _j = cand_het_priority[j];
    //         if (cand_het_pos[_i].count < cand_het_pos[_j].count || (cand_het_pos[_i].count == cand_het_pos[_j].count && cand_het_pos[_i].var_type > cand_het_pos[_j].var_type)) {
    //             printf("Swap: %d (%d, %d) <-> %d (%d, %d)\n", i, cand_het_pos[_i].count, cand_het_pos[_i].var_type, j, cand_het_pos[_j].count, cand_het_pos[_j].var_type);
    //             int tmp = cand_het_priority[i];
    //             cand_het_priority[i] = cand_het_priority[j];
    //             cand_het_priority[j] = tmp;
    //         }
    //     }
    // }
    // print cand_het_pos
    // for (i = 0; i < n_het_pos; ++i) {
    //     int het_i = cand_het_priority[i];
    //     fprintf(stderr, "cand het pos: %d, %d, count: %d, ", i, cand_het_pos[het_i].pos, cand_het_pos[het_i].count);
    //     for (j = 0; j < cand_het_pos[het_i].n_uniq_alles; ++j) {
    //         fprintf(stderr, "%d (%d): ", cand_het_pos[het_i].alle_bases[j], cand_het_pos[het_i].n_clu_reads[cand_het_pos[het_i].alle_bases[j]]);
    //         for (k = 0; k < cand_het_pos[het_i].n_clu_reads[cand_het_pos[het_i].alle_bases[j]]; ++k) {
    //             if (k != 0) fprintf(stderr, ",");
    //             fprintf(stderr, "%d", cand_het_pos[het_i].clu_read_ids[cand_het_pos[het_i].alle_bases[j]][k]);
    //         } fprintf(stderr, "; ");
    //     } fprintf(stderr, "\n");
    // }
    for (j = 0; j < m+1; ++j) free(clu_read_ids[j]);
    free(clu_read_ids); free(n_clu_reads); free(alleles); free(allele_to_idx); free(allele_apperance); free(depth);
    return n_het_pos;
}

int abpoa_collect_msa_dis1(uint8_t **msa, int msa_l, int n_seq, int m, cand_het_pos_t *cand_het_pos, int n_het_pos, int i, int j) {
    int dis = 0;
    for (int k = 0; k < n_het_pos; ++k) {
        int pos = cand_het_pos[k].pos;
    // for (int pos = 0; pos < msa_l; ++pos) {
        int var_weight = 1; // Gap
        if (cand_het_pos[k].var_type == 0) var_weight = 2; // SNP
        uint8_t base1 = msa[i][pos], base2 = msa[j][pos];
        int base1_is_valid = 0, base2_is_valid = 0;
        for (int i = 0; i < cand_het_pos[k].n_uniq_alles; ++i) {
            if (base1 == cand_het_pos[k].alle_bases[i]) base1_is_valid = 1;
            if (base2 == cand_het_pos[k].alle_bases[i]) base2_is_valid = 1;
        }
        if (base1_is_valid == 0 || base2_is_valid == 0) continue;
        if (base1 != base2) dis += var_weight*cand_het_pos[k].count;
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
    // print dis_matrix
    // for (i = 0; i < n_seq; ++i) {
    //     for (j = 0; j < n_seq; ++j) {
    //         printf("%d vs %d: %d\n", i, j, dis_matrix[i][j]);
    //     }
    //     printf("\n");
    // }
    return dis_matrix;
}

// pick two reads from two different clusters with max_dis
int abpoa_collect_2medoids(cand_het_pos_t *cand_het_pos, int het_i, int **dis_matrix, int n_seq, int *med) {
    int max_dis = 0, max_i = -1, max_j = -1;
    for (int i = 0; i < cand_het_pos[het_i].n_uniq_alles-1; ++i) { // i' cluster
        int allele_i = cand_het_pos[het_i].alle_bases[i];
        for (int j = i+1; j < cand_het_pos[het_i].n_uniq_alles; ++j) { // j' cluster
            int allele_j = cand_het_pos[het_i].alle_bases[j];
            // printf("clu %d: %c vs clu %d: %c\n", i, "ACGTN-"[allele_i], j, "ACGTN-"[allele_j]);
            // loop over all reads in i' cluster and j' cluster
            for (int r1 = 0; r1 < cand_het_pos[het_i].n_clu_reads[allele_i]; ++r1) {
                int read_i = cand_het_pos[het_i].clu_read_ids[allele_i][r1];
                for (int r2 = 0; r2 < cand_het_pos[het_i].n_clu_reads[allele_j]; ++r2) {
                    int read_j = cand_het_pos[het_i].clu_read_ids[allele_j][r2];
                    if (dis_matrix[read_i][read_j] > max_dis) {
                        max_dis = dis_matrix[read_i][read_j];
                        max_i = read_i; max_j = read_j;
                    }
                }
            }
        }
    }
    if (max_dis > 0) {
        med[0] = max_i; med[1] = max_j;
        // printf("medoids: %d, %d\n", max_i, max_j);
        return 2;
    } else return 0;
}

int get_partition_index(cand_het_pos_t *cand_het_pos, int het_i, int read_i) {
    int partition_index = 0;
    for (int k = 0; k <= het_i; ++k) {
        partition_index *= (cand_het_pos[k].n_uniq_alles + 1);
        partition_index += cand_het_pos[k].read_id_to_allele_idx[read_i] + 1;
    }
    return partition_index;
}

// 1) if any {0 X ... X het_i}'s clusters is empty, pick one read from the empty cluster with max(min_dis) to existing medoids
// 2) if all {0 X ... X het_i}'s clusters are used, within each cluster, pick one read with max(min_dis) to existing medoids, pick the one with max_dis among the all clusters
int abpoa_collect_1medoid(cand_het_pos_t *cand_het_pos, int het_i, int **dis_matrix, int n_seq, int *med, int n_medoids) {
    assert(n_medoids > 0);
    int n_partitions = 1;
    for (int k = 0; k <= het_i; ++k) {
        n_partitions *= (cand_het_pos[k].n_uniq_alles + 1); // +1 for not covered
    }
    int *partition_counts = (int*)calloc(n_partitions, sizeof(int));
    for (int i = 0; i < n_seq; ++i) {
        int partition_index = get_partition_index(cand_het_pos, het_i, i);
        partition_counts[partition_index] += 1;
    }

    int found = 0;
    // check if any cluster is not used by existing medoids
    int max_dis = 0, max_read_i = -1, max_n_partition_total_count = -1;
    for (int i = 0; i < n_seq; ++i) {
        found = 0;
        int read_partition_index = get_partition_index(cand_het_pos, het_i, i);
        int read_partition_total_count = partition_counts[read_partition_index];
        for (int j = 0; j < n_medoids; ++j) {
            int med_partition_index = get_partition_index(cand_het_pos, het_i, med[j]);
            if (read_partition_index == med_partition_index) {
                found = 1; break;
            }
        }
        
        if (found == 0) { // case 1: read i's partition is not used by any existing medoids
            int min_dis = INT_MAX;
            for (int j = 0; j < n_medoids; ++j) { // min among all medoids
                if (dis_matrix[i][med[j]] < min_dis) min_dis = dis_matrix[i][med[j]];
            }
            if (read_partition_total_count > max_n_partition_total_count || (read_partition_total_count == max_n_partition_total_count && min_dis > max_dis)) { // prefer smaller n_null_allele
                // printf("min: %d, read: %d, n_partition_total_count: %d\n", min_dis, i, read_partition_total_count);
                max_dis = min_dis; max_read_i = i;
                max_n_partition_total_count = read_partition_total_count;
            }
        }
    }
    if (max_read_i == -1) { // case 2: all clusters are used, pick one read from each cluster, then pick the one with max(min_dis)
        for (int i = 0; i < cand_het_pos[het_i].n_uniq_alles; ++i) {
            int allele = cand_het_pos[het_i].alle_bases[i];
            for (int r = 0; r < cand_het_pos[het_i].n_clu_reads[allele]; ++r) {
                int read_i = cand_het_pos[het_i].clu_read_ids[allele][r];
                int min_dis = INT_MAX; int skip = 0;
                for (int j = 0; j < n_medoids; ++j) { // min among all medoids
                    if (med[j] == read_i) {
                        skip = 1;
                        // fprintf(stderr, "Error: medoid %d already in the cluster\n", read_i);
                        continue;
                    }
                    if (dis_matrix[read_i][med[j]] < min_dis) {
                        min_dis = dis_matrix[read_i][med[j]];
                    }
                }
                if (min_dis > max_dis && !skip) {
                    max_dis = min_dis;
                    max_read_i = read_i;
                }
            }
        }
    }
    free(partition_counts);
    if (max_read_i != -1) {
        med[n_medoids] = max_read_i;
        // printf("medoid: %d\n", max_read_i);
        return 1;
    } else return 0;
}

// return up to n_uniq_alleles reads with max dis, each come from one of the cluster defined by het_i
int abpoa_collect_multi_medoids(cand_het_pos_t *cand_het_pos, int het_i, int **dis_matrix, int n_seq, int max_n_cons, int *med, int n_medoids) {
    int n_to_collect = MIN_OF_TWO(cand_het_pos[het_i].n_uniq_alles, max_n_cons);
    while (n_medoids < n_to_collect) {
        int new_medoids = 0;
        if (n_medoids == 0) {
            new_medoids = abpoa_collect_2medoids(cand_het_pos, het_i, dis_matrix, n_seq, med);
        } else {
            new_medoids = abpoa_collect_1medoid(cand_het_pos, het_i, dis_matrix, n_seq, med, n_medoids);
        }
        if (new_medoids == 0) break;
        n_medoids += new_medoids;
    }
    return n_medoids;
}

// allow >2 clusters
int abpoa_init_kmedoids(cand_het_pos_t *cand_het_pos, int *cand_het_priority, int n_het_pos, int **dis_matrix, int n_seq, int max_n_cons, int *medoids) {
    assert(max_n_cons >= 2);
    int n_medoids = 0, het_i = 0;
    while (n_medoids < max_n_cons) {
        if (n_medoids == 0) n_medoids += abpoa_collect_multi_medoids(cand_het_pos, cand_het_priority[het_i], dis_matrix, n_seq, max_n_cons, medoids, n_medoids);
        else n_medoids += abpoa_collect_1medoid(cand_het_pos, cand_het_priority[het_i], dis_matrix, n_seq, medoids, n_medoids);
        het_i++;
        if (het_i >= n_het_pos) break;
    }
    return n_medoids;
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

// XXX swap medoids with non-medoids ???
int abpoa_update_kmedoids(int **dis_matrix, int n_seq, int max_n_cons, int **medoids, int **clu_reads, int *n_clu_seqs, int verbose) {
    int *new_medoids = (int*)malloc(max_n_cons * sizeof(int));
    for (int i = 0; i < max_n_cons; ++i) new_medoids[i] = -1;
    memset(n_clu_seqs, 0, sizeof(int) * max_n_cons);
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
        // if (tied != 0 || min_clu == -1) continue;
        if (min_clu == -1) continue;
        if (tied == 1) {
            if (n_clu_seqs[0] < n_clu_seqs[1]) min_clu = 0;
            else min_clu = 1;
        }
        clu_reads[min_clu][n_clu_seqs[min_clu]++] = i;
    }
    // printf("Clusters after assignment:\n");
    // for (int i = 0; i < max_n_cons; ++i) {
    //     printf("Cluster %d: ", i);
    //     for (int j = 0; j < n_clu_seqs[i]; ++j) {
    //         printf("%d ", clu_reads[i][j]);
    //     }
    //     printf("\n");
    // }
    abpoa_collect_kmedoids0(dis_matrix, n_seq, max_n_cons, clu_reads, n_clu_seqs, new_medoids);
    // print clu_reads and new_medoids
    // printf("Clusters after updating medoids:\n");
    // for (int i = 0; i < max_n_cons; ++i) {
    //     printf("Cluster %d: ", i);
    //     for (int j = 0; j < n_clu_seqs[i]; ++j) {
    //         printf("%d ", clu_reads[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("New Medoids: ");
    // for (int i = 0; i < max_n_cons; ++i) {
    //     printf("%d ", new_medoids[i]);
    // }
    // printf("\n");
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

int abpoa_clu_reads_kmedoids(cand_het_pos_t *cand_het_pos, int *cand_het_priority, int n_het_pos, int **dis_matrix, int n_seq, int min_het, int max_n_cons, uint64_t ***clu_read_ids, int verbose) {
    // iterative update clusters until convergence
    int *medoids = (int*)_err_malloc(sizeof(int) * max_n_cons);
    int **clu_reads = (int**)_err_malloc(sizeof(int*) * max_n_cons);
    int *n_clu_seqs = (int*)_err_malloc(sizeof(int) * max_n_cons);
    for (int i = 0; i < max_n_cons; ++i) {
        clu_reads[i] = (int*)_err_malloc(n_seq * sizeof(int));
    }
    int to_collect_medoids = max_n_cons, n_clusters = 1;
    while (1) {
        if (abpoa_init_kmedoids(cand_het_pos, cand_het_priority, n_het_pos, dis_matrix, n_seq, to_collect_medoids, medoids) <= 0) break; //return 0;
        int iter = 0;
        while (1) {
            int changed = abpoa_update_kmedoids(dis_matrix, n_seq, to_collect_medoids, &medoids, clu_reads, n_clu_seqs, verbose);
            iter++;
            if (changed == 0 || iter >= 10) break;
        }
        int n_clu = 0, n_clustered_reads = 0;
        for (int i = 0; i < to_collect_medoids; ++i) {
            if (n_clu_seqs[i] >= min_het) n_clu++;
            n_clustered_reads += n_clu_seqs[i];
        }
        if (n_clu != to_collect_medoids || n_clustered_reads < ceil(n_seq * 0.8)) {
            to_collect_medoids--;
            if (to_collect_medoids < 2) break;
        } else {
            n_clusters = n_clu;
            break;
        }
    }
    // printf("Final clusters: %d\n", n_clusters);
    if (n_clusters != 1) {
        (*clu_read_ids) = (uint64_t**)_err_malloc(sizeof(uint64_t*) * n_clusters);
        int read_id_n = (n_seq-1)/64+1;
        for (int i = 0; i < n_clusters; ++i) {
            (*clu_read_ids)[i] = (uint64_t*)_err_calloc(read_id_n, sizeof(uint64_t));
            for (int j = 0; j < n_clu_seqs[i]; ++j) {
                int read_i = clu_reads[i][j];
                (*clu_read_ids)[i][read_i/64] |= (1ULL << (read_i & 0x3f));
            }
        }
    }
    for (int i = 0; i < max_n_cons; ++i) free(clu_reads[i]);
    free(clu_reads); free(n_clu_seqs); free(medoids);
    return n_clusters;
}

int abpoa_multip_read_clu_kmedoids(abpoa_graph_t *abg, abpoa_para_t *abpt, int src_id, int sink_id, int n_seq, int m, int max_n_cons, double min_freq, uint64_t ***clu_read_ids, int verbose) {
    abpoa_set_msa_rank(abg, src_id, sink_id);
    int n_clu;
    uint8_t **msa = (uint8_t**)_err_malloc(sizeof(uint8_t*) * n_seq);
    int msa_l = abpoa_collect_msa(abg, abpt, msa, n_seq); // abg->node_id_to_msa_rank[sink_id]-1, 
    // abpoa_collect_allele_freq(msa, msa_l, n_seq, m, verbose);
    int min_w = MAX_OF_TWO(2, ceil(n_seq * min_freq));
    cand_het_pos_t *cand_het_pos = (cand_het_pos_t*)_err_malloc(msa_l * sizeof(cand_het_pos_t));
    int *cand_het_priority = (int*)_err_malloc(msa_l * sizeof(int));
    int n_het_pos = abpoa_collect_cand_het_pos(msa, msa_l, n_seq, m, min_w, cand_het_pos, cand_het_priority, verbose);
    // print cand_het_pos
    // printf("n_het_pos: %d\n", n_het_pos);
    // for (int j = 0; j < n_het_pos; ++j) {
    //     int i = cand_het_priority[j];
    //     fprintf(stdout, "pos: %d, var_type: %d, count: %d, n_uniq_alles: %d, alle_bases: ", cand_het_pos[i].pos, cand_het_pos[i].var_type, cand_het_pos[i].count, cand_het_pos[i].n_uniq_alles);
    //     for (int k = 0; k < cand_het_pos[i].n_uniq_alles; ++k) {
    //         int allele = cand_het_pos[i].alle_bases[k];
    //         fprintf(stdout, "%c: ", "ACGTN-"[allele]);
    //         for (int l = 0; l < cand_het_pos[i].n_clu_reads[allele]; ++l) {
    //             fprintf(stdout, "%d ", cand_het_pos[i].clu_read_ids[allele][l]);
    //         } fprintf(stdout, "; ");
    //     } fprintf(stdout, "\n");
    // }
    // for (int i = 0; i < n_seq; ++i) {
    //     fprintf(stdout, "read %d: ", i);
    //     for (int k = 0; k < 3; ++k) {
    //         int j = cand_het_priority[k];
    //         int allele_idx = cand_het_pos[j].read_id_to_allele_idx[i];
    //         if (allele_idx != -1)
    //             fprintf(stdout, "%c; ", "ACGTN-"[cand_het_pos[j].alle_bases[allele_idx]]);
    //     } fprintf(stdout, "\n");
    // }
    if (n_het_pos < 1) n_clu = 1;
    else {
        int **dis_matrix = abpoa_collect_msa_dis_matrix(msa, msa_l, n_seq, m, cand_het_pos, n_het_pos, verbose);
        n_clu = abpoa_clu_reads_kmedoids(cand_het_pos, cand_het_priority, n_het_pos, dis_matrix, n_seq, min_w, max_n_cons, clu_read_ids, verbose);
        for (int i = 0; i < n_seq; ++i) free(dis_matrix[i]); free(dis_matrix);
    }
    for (int i = 0; i < n_het_pos; ++i) {
        free(cand_het_pos[i].alle_bases);
        for (int j = 0; j < m+1; ++j) free(cand_het_pos[i].clu_read_ids[j]);
        free(cand_het_pos[i].clu_read_ids); free(cand_het_pos[i].n_clu_reads); free(cand_het_pos[i].read_id_to_allele_idx);
    } free(cand_het_pos); free(cand_het_priority);
    for (int i = 0; i < n_seq; ++i) free(msa[i]); free(msa);
    return n_clu;
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
        abpoa_most_frequent(abg, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, n_clu, read_ids_n, clu_read_ids, abc);
    if (n_clu > 1) {
        for (i = 0; i < n_clu; ++i) free(clu_read_ids[i]); free(clu_read_ids);
    }
    abg->is_called_cons = 1; free(out_degree);
}
