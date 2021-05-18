#include <stdio.h>
#include <stdlib.h>
#include "abpoa_align.h"
#include "abpoa_seq.h"
#include "simd_abpoa_align.h"
#include "kdq.h"

char LogTable65536[65536];
char bit_table16[65536];

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
    uint32_t t, tt;
    if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
    return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

void set_65536_table(void) {
    int i;
    for (i = 0; i < 65536; ++i) {
        LogTable65536[i] = ilog2_32(i);
    }
}

void set_bit_table16(void) {
    int i; bit_table16[0] = 0;
    for (i = 0; i != 65536; ++i) bit_table16[i] = (i&1) + bit_table16[i>>1];
}

#define get_bit_cnt4(table, b) (table[(b)&0xffff] + table[(b)>>16&0xffff] + table[(b)>>32&0xffff] + table[(b)>>48&0xffff])

static inline int ilog2_64(uint64_t v) {
    uint64_t t, tt;
    if ((tt = v >> 32)) return (t = tt >> 16) ? 48 + LogTable65536[t] : 32 + LogTable65536[tt];
    return (t = v>>16) ? 16 + LogTable65536[t] : LogTable65536[v];
}

KDQ_INIT(int)
#define kdq_int_t kdq_t(int)

abpoa_node_t *abpoa_init_node(int n) {
    abpoa_node_t *node = (abpoa_node_t*)_err_calloc(n, sizeof(abpoa_node_t));
    return node;
}

void abpoa_set_graph_node(abpoa_graph_t *abg, int node_i) {
    abg->node[node_i].node_id = node_i;
    abg->node[node_i].in_edge_n = 0; abg->node[node_i].in_edge_m = 0;
    abg->node[node_i].out_edge_n = 0; abg->node[node_i].out_edge_m = 0;
    abg->node[node_i].aligned_node_n = 0; abg->node[node_i].aligned_node_m = 0;
    abg->node[node_i].read_ids_n = 0;
}

void abpoa_free_node(abpoa_node_t *node, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (node[i].in_edge_m > 0) free(node[i].in_id);
        if (node[i].out_edge_m > 0) { 
            free(node[i].out_id); free(node[i].out_weight);
            //if (abpt->use_read_ids) {
            if (node[i].read_ids_n > 0)
                free(node[i].read_ids);
            //}
        }
        if (node[i].aligned_node_m > 0) free(node[i].aligned_node_id);
    }
    free(node);
}

// 0: in_edge, 1: out_edge
abpoa_graph_t *abpoa_realloc_graph_edge(abpoa_graph_t *abg, int io, int id) {
    if (io == 0) {
        _uni_realloc(abg->node[id].in_id, abg->node[id].in_edge_n, abg->node[id].in_edge_m, int);
    } else {
        int edge_m = abg->node[id].out_edge_m;
        _uni_realloc(abg->node[id].out_id, abg->node[id].out_edge_n, edge_m, int);
        _uni_realloc(abg->node[id].out_weight, abg->node[id].out_edge_n, abg->node[id].out_edge_m, int);
    }
    return abg;
}

abpoa_graph_t *abpoa_realloc_graph_node(abpoa_graph_t *abg) {
    if (abg->node_m <= 0) {
        abg->node_m = 1;
        abg->node = (abpoa_node_t*)_err_calloc(1, sizeof(abpoa_node_t));
    }
    if (abg->node_n == abg->node_m) {
        int i;
        abg->node_m <<= 1;
        abg->node = (abpoa_node_t*)_err_realloc(abg->node, abg->node_m * sizeof(abpoa_node_t));
        for (i = abg->node_m >> 1; i < abg->node_m; ++i) {
            abpoa_set_graph_node(abg, i);
        }
    }
    return abg;
}

abpoa_graph_t *abpoa_init_graph(void) {
    abpoa_graph_t *abg = (abpoa_graph_t*)_err_malloc(sizeof(abpoa_graph_t));
    abg->node_n = 2, abg->node_m = 2, abg->index_rank_m = 0;
    abg->node = abpoa_init_node(2);
    abg->node[0].node_id = 0; abg->node[1].node_id = 1;
    abg->node[0].read_ids_n = 0; abg->node[1].read_ids_n = 0;
    abg->is_topological_sorted = abg->is_called_cons = 0;
    abg->node_id_to_index = NULL; abg->index_to_node_id = NULL; abg->node_id_to_msa_rank = NULL;
    abg->node_id_to_max_pos_left = NULL; abg->node_id_to_max_pos_right = NULL; abg->node_id_to_max_remain = NULL;
    return abg;
}

void abpoa_free_graph(abpoa_graph_t *abg) {
    if (abg->node_m > 0) abpoa_free_node(abg->node, abg->node_m);

    if (abg->node_n > 0) {
        free(abg->index_to_node_id);
        free(abg->node_id_to_index);
        if (abg->node_id_to_msa_rank) free(abg->node_id_to_msa_rank);

        if (abg->node_id_to_max_pos_left) free(abg->node_id_to_max_pos_left);
        if (abg->node_id_to_max_pos_right) free(abg->node_id_to_max_pos_right);
        if (abg->node_id_to_max_remain) free(abg->node_id_to_max_remain);
    }
    free(abg);
}

abpoa_t *abpoa_init(void) {
    abpoa_t *ab = (abpoa_t*)_err_malloc(sizeof(abpoa_t));
    ab->abg = abpoa_init_graph();
    ab->abs = abpoa_init_seq();
    ab->abm = abpoa_init_simd_matrix();
    return ab;
}

void abpoa_free(abpoa_t *ab) {
    abpoa_free_graph(ab->abg);
    abpoa_free_seq(ab->abs);
    abpoa_free_simd_matrix(ab->abm);
    free(ab);
}

void abpoa_BFS_set_node_index(abpoa_graph_t *abg, int src_id, int sink_id) {
    int *id, cur_id, out_id, aligned_id;
    int index = 0, q_size, new_q_size;

    int *in_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    int i, j;
    for (i = 0; i < abg->node_n; ++i) in_degree[i] = abg->node[i].in_edge_n;

    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, src_id); q_size = 1; new_q_size = 0; // node[q.id].in_degree equals 0
    while (q_size > 0) {
        if ((id = kdq_shift_int(q)) == 0) err_fatal_simple("Error in queue.");
        cur_id = *id;
        abg->index_to_node_id[index] = cur_id;
        abg->node_id_to_index[cur_id] = index++;

        if (cur_id == sink_id) {
            kdq_destroy_int(q); free(in_degree);
            return;
        }
        for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
            out_id = abg->node[cur_id].out_id[i];
            if (--in_degree[out_id] == 0) {
                for (j = 0; j < abg->node[out_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[out_id].aligned_node_id[j];
                    if (in_degree[aligned_id] != 0) goto next_out_node;
                }
                kdq_push_int(q, out_id);
                ++new_q_size;
                for (j = 0; j < abg->node[out_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[out_id].aligned_node_id[j];
                    kdq_push_int(q, aligned_id);
                    ++new_q_size;
                }
            }
next_out_node:;
        }
        if (--q_size == 0) {
            q_size = new_q_size;
            new_q_size = 0;
        }
    }
    err_fatal_simple("Failed to set node index.");
}

void abpoa_BFS_set_node_remain(abpoa_graph_t *abg, int src_id, int sink_id) {
    int *id, cur_id, i, out_id, in_id;

    int *out_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    for (i = 0; i < abg->node_n; ++i) {
        out_degree[i] = abg->node[i].out_edge_n;
        abg->node_id_to_max_remain[i] = 0;
    }

    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, sink_id); // node[q.id].in_degree equals 0
    abg->node_id_to_max_remain[sink_id] = -1; // XXX not 0
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;

        // all out_id of cur_id have beed visited
        // max weight out_id
        if (cur_id != sink_id) {
            int max_w=-1, max_id=-1;
            for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                out_id = abg->node[cur_id].out_id[i];
                if (abg->node[cur_id].out_weight[i] > max_w) {
                    max_w = abg->node[cur_id].out_weight[i];
                    max_id = out_id;
                }
            }
            abg->node_id_to_max_remain[cur_id] = abg->node_id_to_max_remain[max_id] + 1;
            // fprintf(stderr, "%d -> %d\n", abg->node_id_to_index[cur_id], abg->node_id_to_max_remain[cur_id]);
        }
        if (cur_id == src_id) {
            kdq_destroy_int(q); free(out_degree);
            return;
        }
        for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
            in_id = abg->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) kdq_push_int(q, in_id);
        }
    }
    err_fatal_simple("Failed to set node remain.");
}

// 1. index_to_node_id
// 2. node_id_to_index
// 3. node_id_to_rank
void abpoa_topological_sort(abpoa_graph_t *abg, abpoa_para_t *abpt) {
    if (abg->node_n <= 0) {
        err_func_format_printf(__func__, "Empty graph.\n");
        return;
    }
    int node_n = abg->node_n;
    if (node_n > abg->index_rank_m) {
        abg->index_rank_m = node_n; kroundup32(abg->index_rank_m);
        // fprintf(stderr, "node_n: %d, index_rank_m: %d\n", node_n, abg->index_rank_m);
        abg->index_to_node_id = (int*)_err_realloc(abg->index_to_node_id, abg->index_rank_m * sizeof(int));
        abg->node_id_to_index = (int*)_err_realloc(abg->node_id_to_index, abg->index_rank_m * sizeof(int));
        if (abpt->out_msa || abpt->cons_agrm == ABPOA_HC || abpt->is_diploid) 
            abg->node_id_to_msa_rank = (int*)_err_realloc(abg->node_id_to_msa_rank, abg->index_rank_m * sizeof(int));
        if (abpt->wb >= 0) {
            abg->node_id_to_max_pos_left = (int*)_err_realloc(abg->node_id_to_max_pos_left, abg->index_rank_m * sizeof(int));
            abg->node_id_to_max_pos_right = (int*)_err_realloc(abg->node_id_to_max_pos_right, abg->index_rank_m * sizeof(int));
            abg->node_id_to_max_remain = (int*)_err_realloc(abg->node_id_to_max_remain, abg->index_rank_m * sizeof(int));
        } else if (abpt->zdrop > 0) {
            abg->node_id_to_max_remain = (int*)_err_realloc(abg->node_id_to_max_remain, abg->index_rank_m * sizeof(int));
        }
    }
    // start from ABPOA_SRC_NODE_ID to ABPOA_SINK_NODE_ID
    abpoa_BFS_set_node_index(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    // init min/max rank
    if (abpt->wb >= 0) {
        int i;
        for (i = 0; i < node_n; ++i) {
            abg->node_id_to_max_pos_right[i] = 0;
            abg->node_id_to_max_pos_left[i] = node_n;
        }
        abpoa_BFS_set_node_remain(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    } else if (abpt->zdrop > 0)
        abpoa_BFS_set_node_remain(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    abg->is_topological_sorted = 1;
}

void abpoa_DFS_set_msa_rank(abpoa_graph_t *abg, int src_id, int sink_id, int *in_degree) {
    // fprintf(stderr, "node_n: %d, m: %d\n", abg->node_n, abg->index_rank_m);
    if (abg->node_n > abg->index_rank_m) {
        int m = abg->node_n; kroundup32(m);
        abg->node_id_to_msa_rank = (int*)_err_realloc(abg->node_id_to_msa_rank, m * sizeof(int));
    }
    int *id, cur_id, i, j, out_id, aligned_id;
    int msa_rank = 0;
    kdq_int_t *q = kdq_init_int();

    // Depth-First-Search
    kdq_push_int(q, src_id); // node[q.id].in_degree equals 0
    abg->node_id_to_msa_rank[src_id] = -1;
    // printf("tot_node_n: %d, node_m: %d\n", abg->node_n, abg->node_m);

    while((id = kdq_pop_int(q)) != 0) {
        cur_id = *id;
        if (abg->node_id_to_msa_rank[cur_id] < 0) {
            abg->node_id_to_msa_rank[cur_id] = msa_rank;
            for (i = 0; i < abg->node[cur_id].aligned_node_n; ++i) {
                aligned_id = abg->node[cur_id].aligned_node_id[i];
                abg->node_id_to_msa_rank[aligned_id] = msa_rank;
            }
            msa_rank++;
        }

        if (cur_id == sink_id) {
            kdq_destroy_int(q);
            abg->is_set_msa_rank = 1;
            return;
        }
        for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
            out_id = abg->node[cur_id].out_id[i];
            if (--in_degree[out_id] == 0) {
                for (j = 0; j < abg->node[out_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[out_id].aligned_node_id[j];
                    if (in_degree[aligned_id] != 0) goto next_out_node;
                }
                kdq_push_int(q, out_id);
                abg->node_id_to_msa_rank[out_id] = -1;
                for (j = 0; j < abg->node[out_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[out_id].aligned_node_id[j];
                    kdq_push_int(q, aligned_id);
                    // printf("aln_id: %d\n", aligned_id);
                    abg->node_id_to_msa_rank[aligned_id] = -1;
                }
            }
next_out_node:;
        }
    }
    err_fatal_simple("Error in set_msa_rank.\n");
}

void abpoa_set_msa_rank(abpoa_graph_t *abg, int src_id, int sink_id) {
    if (abg->is_set_msa_rank == 0) {
        int i, *in_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
        for (i = 0; i < abg->node_n; ++i) in_degree[i] = abg->node[i].in_edge_n;
        abpoa_DFS_set_msa_rank(abg, src_id, sink_id, in_degree);
        free(in_degree);
    }
}

void abpoa_set_row_column_weight(abpoa_graph_t *abg, int **rc_weight, int **msa_node_id) {
    int i, k, rank, aligned_id;
    uint64_t b;
    for (i = 2; i < abg->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(abg, i);
        for (k = 0; k < abg->node[i].aligned_node_n; ++k) {
            aligned_id = abg->node[i].aligned_node_id[k];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
        }
        // assign seq
        for (k = 0; k < abg->node[i].read_ids_n; ++k) {
            b = abg->node[i].read_ids[k];
            rc_weight[rank-1][abg->node[i].base] += get_bit_cnt4(bit_table16, b);
        }
        rc_weight[rank-1][4] -= rc_weight[rank-1][abg->node[i].base];
        msa_node_id[rank-1][abg->node[i].base] = i;
    }
}

void abpoa_set_row_column_ids_weight(abpoa_graph_t *abg, uint64_t ***read_ids, int **rc_weight, int **msa_node_id, int msa_l, int n_seq, int read_ids_n) {
    int i, j, k, rank, aligned_id;
    uint64_t b, one = 1, *whole_read_ids = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    for (i = 0; i < n_seq; ++i) {
        j = i / 64; b = i & 0x3f;
        whole_read_ids[j] |= (one << b);
    }
    for (i = 0; i < msa_l; ++i) {
        for (j = 0; j < read_ids_n; ++j) {
            read_ids[i][4][j] = whole_read_ids[j];
        }
    }
    for (i = 2; i < abg->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(abg, i);
        for (k = 0; k < abg->node[i].aligned_node_n; ++k) {
            aligned_id = abg->node[i].aligned_node_id[k];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
        }
        // assign seq
        for (k = 0; k < abg->node[i].read_ids_n; ++k) {
            b = abg->node[i].read_ids[k];
            rc_weight[rank-1][abg->node[i].base] += get_bit_cnt4(bit_table16, b);
            read_ids[rank-1][abg->node[i].base][k] = b;
            read_ids[rank-1][4][k] ^= b;
        }
        rc_weight[rank-1][4] -= rc_weight[rank-1][abg->node[i].base];
        msa_node_id[rank-1][abg->node[i].base] = i;
    }
    free(whole_read_ids);
}

void abpoa_heaviest_column_multip_consensus(uint64_t ***read_ids, int **cluster_ids, int *cluster_ids_n, int multip, int read_ids_n, int msa_l, FILE *out_fp, uint8_t ***_cons_seq, int **_cons_l, int *_cons_n) {
    int i, j, k, m, w, cnt, max_w, max_base, gap_w;
    uint64_t *read_ids_mask = (uint64_t*)_err_malloc(read_ids_n * sizeof(uint64_t)), one=1, b;
    uint8_t *cons_seq = (uint8_t*)_err_malloc(sizeof(uint8_t) * msa_l); int cons_l;
    int read_id, n;

    if (_cons_n) {
        *_cons_n = multip; 
        (*_cons_l) = (int*)_err_malloc(sizeof(int) * multip);
        (*_cons_seq) = (uint8_t **)_err_malloc(sizeof(uint8_t*) * multip);
    }
    for (i = 0; i < multip; ++i) {
        cons_l = 0;
        if (out_fp) fprintf(out_fp, ">Consensus_sequence_%d_%d\n", i+1, cluster_ids_n[i]);
        for (j = 0; j < read_ids_n; ++j) read_ids_mask[j] = 0;
        for (j = 0; j < cluster_ids_n[i]; ++j) {
            read_id = cluster_ids[i][j];
            n = read_id / 64;
            read_ids_mask[n] |= (one << (read_id & 0x3f));
        }
        for (j = 0; j < msa_l; ++j) {
            max_w = 0, max_base = 5, gap_w = cluster_ids_n[i]; 
            for (k = 0; k < 4; ++k) {
                w = 0;
                for (m = 0; m < read_ids_n; ++m) {
                    b = read_ids[j][k][m] & read_ids_mask[m];
                    cnt = get_bit_cnt4(bit_table16, b);
                    w += cnt, gap_w -= cnt;
                }
                if (w > max_w) {
                    max_w = w; max_base = k;
                }
            }
            if (max_w > gap_w) {
                cons_seq[cons_l++] = max_base;
            }
        }
        if (out_fp) {
            for (j = 0; j < cons_l; ++j) fprintf(out_fp, "%c", "ACGT"[cons_seq[j]]); fprintf(out_fp, "\n");
        }
        if (_cons_n) {
            (*_cons_l)[i] = cons_l;
            (*_cons_seq)[i] = (uint8_t*)_err_malloc(sizeof(uint8_t) * cons_l);
            for (j = 0; j < cons_l; ++j) (*_cons_seq)[i][j] = cons_seq[j];
        }
    }
    free(read_ids_mask); free(cons_seq);
}

int output_consensus(abpoa_graph_t *abg, int src_id, int sink_id, FILE *out_fp) {
    int cons_l = 0;
    fprintf(out_fp, ">Consensus_sequence\n");
    int id = abg->node[src_id].max_out_id;
    while (id != sink_id) {
        fprintf(out_fp, "%c", "ACGTN"[abg->node[id].base]);
        id = abg->node[id].max_out_id;
        cons_l++;
    } fprintf(out_fp, "\n");
    return cons_l;
}

int abpoa_store_consensus(abpoa_graph_t *abg, int src_id, int sink_id, uint8_t ***cons_seq, int **cons_l) {
    *cons_seq = (uint8_t**)_err_malloc(sizeof(uint8_t*));
    *cons_l = (int*)_err_malloc(sizeof(int));
    (*cons_seq)[0] = (uint8_t*)_err_malloc(sizeof(uint8_t) * abg->node_n);
    int id = abg->node[src_id].max_out_id, i = 0;
    while (id != sink_id) {
        (*cons_seq)[0][i++] = abg->node[id].base;
        id = abg->node[id].max_out_id;
    }
    (*cons_l)[0] = i;
    return i;
}

void abpoa_generate_consensus_cov(abpoa_graph_t *abg, int src_id, int sink_id, int ***cons_cov) {
    *cons_cov = (int**)_err_malloc(sizeof(int*));
    (*cons_cov)[0] = (int*)_err_malloc(sizeof(int) * abg->node_n);
    int id = abg->node[src_id].max_out_id;
    int left_w, right_w;
    int i, j, in_id, cons_i = 0;
    while (id != sink_id) {
        // for each id: get max{left_weigth, right_weight}
        left_w = right_w = 0;
        for (i = 0; i < abg->node[id].in_edge_n; ++i) {
            in_id = abg->node[id].in_id[i];
            for (j = 0; j < abg->node[in_id].out_edge_n; ++j) {
                if (abg->node[in_id].out_id[j] == id) {
                    left_w += abg->node[in_id].out_weight[j];
                }
            }
        }
        for (i = 0; i < abg->node[id].out_edge_n; ++i) {
            right_w += abg->node[id].out_weight[i];
        }
        (*cons_cov)[0][cons_i++] = MAX_OF_TWO(left_w, right_w);
        id = abg->node[id].max_out_id;
    }
}

void abpoa_traverse_min_flow(abpoa_graph_t *abg, int src_id, int sink_id, int *out_degree) {
    int *id, i, cur_id, in_id, out_id, max_id; int max_w, max_out_i, min_w;
    kdq_int_t *q = kdq_init_int(); kdq_push_int(q, sink_id);
    int *max_out_weight = (int*)_err_malloc(abg->node_n * sizeof(int));
    // reverse Breadth-First-Search
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == sink_id) {
            abg->node[sink_id].max_out_id = -1;
            max_out_weight[sink_id] = INT32_MAX;
        } else {
            max_w = INT32_MIN, max_id = -1, max_out_i = -1, min_w = INT32_MAX;
            for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                out_id = abg->node[cur_id].out_id[i];
                min_w = MIN_OF_TWO(max_out_weight[out_id], abg->node[cur_id].out_weight[i]);
                if (max_w < min_w) {
                    max_w = min_w;
                    max_id = out_id;
                    max_out_i = i;
                } else if (max_w == min_w && abg->node[cur_id].out_weight[max_out_i] <= abg->node[cur_id].out_weight[i]) {
                    max_id = out_id;
                    max_out_i = i;
                }
            }
            abg->node[cur_id].max_out_id = max_id;
            max_out_weight[cur_id] = max_w;
        }

        if (cur_id == src_id) { 
            kdq_destroy_int(q); 
            goto MF_CONS; 
        }
        for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
            in_id = abg->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) kdq_push_int(q, in_id);
        }
    }
MF_CONS:
    free(max_out_weight);
}

// heaviest_bundling
// 1. argmax{cur->weight}
// 2. argmax{out_node->weight}
void abpoa_heaviest_bundling(abpoa_graph_t *abg, int src_id, int sink_id, int *out_degree, int ***cons_cov) {
    int *id, i, cur_id, in_id, out_id, max_id; int max_w, out_w;

    int *score = (int*)_err_malloc(abg->node_n * sizeof(int));
    kdq_int_t *q = kdq_init_int();
    kdq_push_int(q, sink_id);
    // reverse Breadth-First-Search
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == sink_id) {
            abg->node[sink_id].max_out_id = -1;
            score[sink_id] = 0;
        } else {
            max_id = -1;
            if (cur_id == src_id) {
                int path_score = -1, path_max_w = -1;
                for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                    out_id = abg->node[cur_id].out_id[i];
                    out_w = abg->node[cur_id].out_weight[i];
                    if (out_w > path_max_w || (out_w == path_max_w && score[out_id] > path_score)) {
                        max_id = out_id;
                        path_score = score[out_id];
                        path_max_w = out_w;
                    }
                }
                abg->node[src_id].max_out_id = max_id;
                kdq_destroy_int(q);
                goto HB_CONS;
            } else {
                max_w = INT32_MIN;
                for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                    out_id = abg->node[cur_id].out_id[i];
                    out_w = abg->node[cur_id].out_weight[i];
                    if (max_w < out_w) {
                        max_w = out_w;
                        max_id = out_id;
                    } else if (max_w == out_w && score[max_id] <= score[out_id]) {
                        max_id = out_id;
                    }
                }
                score[cur_id] = max_w + score[max_id];
                abg->node[cur_id].max_out_id = max_id;
            }
        }
        for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
            in_id = abg->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) 
                kdq_push_int(q, in_id);
        }
    }
HB_CONS:
    if (cons_cov != NULL) abpoa_generate_consensus_cov(abg, src_id, sink_id, cons_cov); 
    free(score);
}

void abpoa_heaviest_column_consensus(abpoa_graph_t *abg, int **rc_weight, int **msa_node_id, int src_id, int sink_id, int msa_l, int n_seq) {
    int i, j, w, max_base, max_w, gap_w;
    int last_id = src_id, cur_id;
    for (i = 0; i < msa_l; ++i) {
        max_w = 0, max_base = 5, gap_w = n_seq;
        for (j = 0; j < 4; ++j) {
            w = rc_weight[i][j];
            if (w > max_w) {
                max_base = j, max_w = w;
            }
            gap_w -= w;
        }
        if (max_w >= gap_w) {
            cur_id = msa_node_id[i][max_base];
            abg->node[last_id].max_out_id = cur_id;
            last_id = cur_id;
        }
    }
    abg->node[last_id].max_out_id = sink_id;
}

void abpoa_heaviest_column(abpoa_graph_t *abg, int src_id, int sink_id, int n_seq) {
    abpoa_set_msa_rank(abg, src_id, sink_id);

    int i, msa_l = abg->node_id_to_msa_rank[sink_id] - 1;
    int **rc_weight = (int**)_err_malloc(sizeof(int*) * msa_l);
    int **msa_node_id = (int**)_err_malloc(sizeof(int*) * msa_l);
    for (i = 0; i < msa_l; ++i) {
        rc_weight[i] = (int*)_err_calloc(5, sizeof(int)); // ACGT
        rc_weight[i][4] = n_seq;
        msa_node_id[i] = (int*)_err_calloc(5, sizeof(int));
    } 
    abpoa_set_row_column_weight(abg, rc_weight, msa_node_id);

    abpoa_heaviest_column_consensus(abg, rc_weight, msa_node_id, src_id, sink_id, msa_l, n_seq);
    for (i = 0; i < msa_l; ++i) { 
        free(rc_weight[i]); free(msa_node_id[i]);
    } free(rc_weight); free(msa_node_id);
}

void add_het_read_ids(int *init, uint64_t **het_read_ids, uint8_t **het_cons_base, uint64_t **read_ids, int het_n, int *multip_i, int read_ids_n) {
    int i, j;
    if (*init) {
        for (i = 0; i < 2; ++i) {
            for (j = 0; j < read_ids_n; ++j) het_read_ids[i][j] = read_ids[multip_i[i]][j];
            het_cons_base[het_n][0] = multip_i[0];
            het_cons_base[het_n][1] = multip_i[1];
        }
        *init = 0;
    } else {
        int ovlp_n, max_ovlp_n, max_clu_i;
        uint64_t *ids, b;
        max_ovlp_n = max_clu_i = -1;
        ids = read_ids[multip_i[0]];
        for (i = 0; i < 2; ++i) { // het_read_ids
            ovlp_n = 0; 
            for (j = 0; j < read_ids_n; ++j) {
                b = het_read_ids[i][j] & ids[j];
                ovlp_n += get_bit_cnt4(bit_table16, b);
            }
            if (ovlp_n > max_ovlp_n) {
                max_ovlp_n = ovlp_n;
                max_clu_i = i;
            }
        }
        // 0 & i, 1 & 1-i
        for (i = 0; i < read_ids_n; ++i) {
            het_read_ids[max_clu_i][i] |= read_ids[multip_i[0]][i];
            het_read_ids[1-max_clu_i][i] |= read_ids[multip_i[1]][i];
        }
        het_cons_base[het_n][max_clu_i] = multip_i[0];
        het_cons_base[het_n][1-max_clu_i] = multip_i[1];
    }
}

int set_clu_read_ids(int **clu_read_ids, int *clu_read_ids_n, int n_seq, double min_freq, uint64_t ***read_ids, uint8_t **het_cons_base, int *het_pos, int het_n) {
    int i, j, seq_i, pos, clu_base;
    int **diff = (int**)_err_malloc(sizeof(int*) * n_seq);
    for (i = 0; i < n_seq; ++i) {
        diff[i] = (int*)_err_malloc(2 *  sizeof(int));
        diff[i][0] = het_n, diff[i][1] = het_n;
    }
    uint64_t seq_b, one = 1; int seq_bit;
    for (seq_i = 0; seq_i < n_seq; ++seq_i) {
        seq_bit = seq_i / 64; seq_b = one << (seq_i & 0x3f); 
        for (i = 0; i < het_n; ++i) {
            pos = het_pos[i];
            for (j = 0; j < 2; ++j) {
                clu_base = het_cons_base[i][j];
                if (seq_b & read_ids[pos][clu_base][seq_bit]) {
                    --diff[seq_i][j]; break;
                }
            }
        }
    }
    for (i = 0; i < n_seq; ++i) {
        if (diff[i][0] < diff[i][1]) {
            clu_read_ids[0][clu_read_ids_n[0]++] = i;
        } else if (diff[i][1] < diff[i][0]) {
            clu_read_ids[1][clu_read_ids_n[1]++] = i;
        }
    }
    for (i = 0; i < n_seq; ++i) free(diff[i]); free(diff);

    int min_n = MIN_OF_TWO(clu_read_ids_n[0], clu_read_ids_n[1]);
    int clu_n = (min_n >= (int)(n_seq * min_freq)) ? 2 : 1; 
    return clu_n;
}

int abpoa_diploid_ids(uint64_t ***read_ids, int **rc_weight, int msa_l, int n_seq, double min_freq, int read_ids_n, int **clu_read_ids, int *clu_read_ids_n) {
    int i, j;
    uint64_t **het_read_ids = (uint64_t**)_err_malloc(sizeof(uint64_t*) * 2 );
    for (i = 0; i < 2; ++i) het_read_ids[i] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    int *het_pos = (int*)_err_malloc(sizeof(int) * msa_l), het_n = 0;
    uint8_t **het_cons_base = (uint8_t**)_err_malloc(sizeof(uint8_t*) * msa_l);
    for (i = 0; i < msa_l; ++i) het_cons_base[i] = (uint8_t*)_err_malloc(sizeof(uint8_t) * 2);
    int min_w = MAX_OF_TWO(1, n_seq * min_freq);
    int init = 1, tmp, multip_i[2], multip_n, w;
    // collect het nodes
    for (i = 0; i < msa_l; ++i) {
        multip_n = 0;
        for (j = 0; j < 5; ++j) {
            w = rc_weight[i][j];
            if (w >= min_w) {
                if (multip_n == 2) {
                    multip_n = 0;
                    break;
                }
                multip_i[multip_n++] = j;
            }
        }
        if (multip_n == 2) {
            if (rc_weight[i][multip_i[0]] < rc_weight[i][multip_i[1]]) {
                tmp = multip_i[0]; multip_i[0] = multip_i[1]; multip_i[1] = tmp;
            }
            // iteratively update read_ids and cons-base for each cluster
            //     read_ids |=
            //     cons_base1[i++] = ; cons_base2[i++] = ;
            add_het_read_ids(&init, het_read_ids, het_cons_base, read_ids[i], het_n, multip_i, read_ids_n);
            het_pos[het_n++] = i;
        }
    }
    int clu_n;
    if (het_n == 0) clu_n = 1;
    else {
        // for each read, determine cluster based on diff with two cons-bases
        // return two group of read_ids
        clu_n = set_clu_read_ids(clu_read_ids, clu_read_ids_n, n_seq, min_freq, read_ids, het_cons_base, het_pos, het_n);
    }
    
    for (i = 0; i < msa_l; ++i) free(het_cons_base[i]);
    for (i = 0; i < 2; ++i) free(het_read_ids[i]); 
    free(het_cons_base); free(het_read_ids); free(het_pos);
    return clu_n;
}

void abpoa_diploid_heaviest_column(abpoa_graph_t *abg, int src_id, int sink_id, int n_seq, double min_freq, FILE *out_fp, uint8_t ***cons_seq, int **cons_l, int *cons_n) {
    abpoa_set_msa_rank(abg, src_id, sink_id);
    int i, j, msa_l = abg->node_id_to_msa_rank[sink_id] - 1;
    int read_ids_n = (n_seq-1)/64+1;
    uint64_t ***read_ids = (uint64_t***)_err_malloc(sizeof(uint64_t**) * msa_l);
    for (i = 0; i < msa_l; ++i) {
        read_ids[i] = (uint64_t**)_err_malloc(sizeof(uint64_t*) * 5);
        for (j = 0; j < 5; ++j) read_ids[i][j] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    }
    int **rc_weight = (int**)_err_malloc(msa_l * sizeof(int*));
    int **msa_node_id = (int**)_err_malloc(msa_l * sizeof(int*));
    for (i = 0; i < msa_l; ++i) {
        rc_weight[i] = (int*)_err_calloc(5, sizeof(int)); // ACGT
        rc_weight[i][4] = n_seq;
        msa_node_id[i] = (int*)_err_calloc(5, sizeof(int)); // ACGT
    } 
    abpoa_set_row_column_ids_weight(abg, read_ids, rc_weight, msa_node_id, msa_l, n_seq, read_ids_n);

    int **clu_read_ids = (int**)_err_malloc(sizeof(int*) * 2), *clu_read_ids_n;
    for (i = 0; i < 2; ++i) clu_read_ids[i] = (int*)_err_malloc(sizeof(int) * n_seq);
    clu_read_ids_n = (int*)_err_calloc(2, sizeof(int));
    int clu_n = abpoa_diploid_ids(read_ids, rc_weight, msa_l, n_seq, min_freq, read_ids_n, clu_read_ids, clu_read_ids_n);
    if (clu_n == 1) {
        abpoa_heaviest_column_consensus(abg, rc_weight, msa_node_id, src_id, sink_id, msa_l, n_seq);
        if (out_fp) output_consensus(abg, src_id, sink_id, out_fp);
        if (cons_n) {
            *cons_n = 1; abpoa_store_consensus(abg, src_id, sink_id, cons_seq, cons_l);
        }
    } else abpoa_heaviest_column_multip_consensus(read_ids, clu_read_ids, clu_read_ids_n, clu_n, read_ids_n, msa_l, out_fp, cons_seq, cons_l, cons_n);

    for (i = 0; i < msa_l; ++i) {
        for (j = 0; j < 5; ++j) {
            free(read_ids[i][j]);
        } 
        free(read_ids[i]); free(rc_weight[i]); free(msa_node_id[i]);
    } free(read_ids); free(rc_weight); free(msa_node_id);
    for (i = 0; i < 2; ++i) free(clu_read_ids[i]); free(clu_read_ids); free(clu_read_ids_n);
}

// should always topological sort first, then generate consensus
int abpoa_generate_consensus(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp, uint8_t ***cons_seq, int ***cons_cov, int **cons_l, int *cons_n) {
    abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) return 0;
    int i, _cons_l = 0, *out_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    int n_seq = ab->abs->n_seq;
    for (i = 0; i < abg->node_n; ++i) {
        out_degree[i] = abg->node[i].out_edge_n;
    }

    if (abpt->is_diploid) {
        abpoa_diploid_heaviest_column(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, n_seq, abpt->min_freq, out_fp, cons_seq, cons_l, cons_n);
    } else {
        if (abpt->cons_agrm == ABPOA_HB) abpoa_heaviest_bundling(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree, cons_cov);
        else if (abpt->cons_agrm == ABPOA_HC) abpoa_heaviest_column(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, n_seq);
        else if (abpt->cons_agrm == ABPOA_MF) abpoa_traverse_min_flow(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree); 
        else err_fatal(__func__, "Unknown consensus calling algorithm: %d.", abpt->cons_agrm);
        if (out_fp) _cons_l = output_consensus(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_fp);
        if (cons_n) {
            *cons_n = 1; 
            _cons_l = abpoa_store_consensus(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, cons_seq, cons_l);
        }
    } 
    abg->is_called_cons = 1; free(out_degree);
    return _cons_l;
}

void abpoa_set_msa_seq(abpoa_node_t node, int rank, uint8_t **msa_seq) {
    int i, b, read_id; uint8_t base = node.base;
    uint64_t num, tmp;
    b = 0;
    for (i = 0; i < node.read_ids_n; ++i) {
        num = node.read_ids[i];
        while (num) {
            tmp = num & -num;
            read_id = ilog2_64(tmp);
            // fprintf(stderr, "%d -> %d\n", node.node_id, read_id);
            msa_seq[b+read_id][rank-1] = base;
            num ^= tmp;
        }
        b += 64;
    }
}

void abpoa_generate_rc_msa(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp, uint8_t ***msa_seq, int *msa_l) {
    abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) return;
    abpoa_set_msa_rank(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);

    abpoa_seq_t *abs = ab->abs;
    int i, j, k, aligned_id, rank, n_seq = abs->n_seq;
    uint8_t **_msa_seq = (uint8_t**)_err_malloc(n_seq * sizeof(uint8_t*));
    int _msa_l = abg->node_id_to_msa_rank[ABPOA_SINK_NODE_ID]-1;

    for (i = 0; i < n_seq; ++i) {
        _msa_seq[i] = (uint8_t*)_err_malloc(_msa_l * sizeof(uint8_t));
        for (j = 0; j < _msa_l; ++j) 
            _msa_seq[i][j] = 5;
    }

    if (out_fp && abpt->out_msa_header == 0) fprintf(out_fp, ">Multiple_sequence_alignment\n");
    for (i = 2; i < abg->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(abg, i);
        for (k = 0; k < abg->node[i].aligned_node_n; ++k) {
            aligned_id = abg->node[i].aligned_node_id[k];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
        }
        // assign seq
        abpoa_set_msa_seq(abg->node[i], rank, _msa_seq);
    }
    if (out_fp) {
        for (i = 0; i < n_seq; ++i) {
            if (abpt->out_msa_header && abs->name[i].l > 0) {
                if (abs->is_rc[i]) fprintf(out_fp, ">%s_reverse_complement\n", abs->name[i].s);
                else fprintf(out_fp, ">%s\n", abs->name[i].s);
            }
            for (j = 0; j < _msa_l; ++j) fprintf(out_fp, "%c", "ACGTN-"[_msa_seq[i][j]]);
            fprintf(out_fp, "\n");
        }
        if (abpt->out_cons) { // RC-MSA for consensus sequence
            if (abpt->out_msa_header && abs->name[0].l > 0) fprintf(out_fp, ">Consensus_sequence\n");

            i = abg->node[ABPOA_SRC_NODE_ID].max_out_id;
            int last_rank = 1;
            while (i != ABPOA_SINK_NODE_ID) {
                rank = abpoa_graph_node_id_to_msa_rank(abg, i);
                for (k = 0; k < abg->node[i].aligned_node_n; ++k) {
                    aligned_id = abg->node[i].aligned_node_id[k];
                    rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
                }
                // last_rank -> rank : -
                for (k = last_rank; k < rank; ++k) fprintf(out_fp, "-");
                // rank : base
                fprintf(out_fp, "%c", "ACGTN"[abg->node[i].base]);
                last_rank = rank+1;
                i = abg->node[i].max_out_id;
            }
            // last_rank -> msa_l:-
            for (k = last_rank; k <= _msa_l; ++k) fprintf(out_fp, "-");
            fprintf(out_fp, "\n");
        }
    }
    if (msa_l) {
        *msa_l = _msa_l; *msa_seq = _msa_seq;
    } else {
        for (i = 0; i < n_seq; ++i) free(_msa_seq[i]); free(_msa_seq);
    }
}

void abpoa_generate_gfa(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp) {
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
                fprintf(out_fp, "S\t%d\t%c\n", cur_id-1, "ACGTN"[abg->node[cur_id].base]);
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
                    num = abg->node[cur_id].read_ids[i];
                    while (num) {
                        tmp = num & -num;
                        read_id = ilog2_64(tmp);
                        read_paths[b+read_id][read_path_i[b+read_id]++] = cur_id-1;
                        num ^= tmp;
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
        abpoa_generate_consensus(ab, abpt, NULL, NULL, NULL, NULL, NULL);
        int id = abg->node[ABPOA_SRC_NODE_ID].max_out_id;
        fprintf(out_fp, "P\tConsensus_sequence\t");
        while (1) {
            fprintf(out_fp, "%d+", id-1);
            id = abg->node[id].max_out_id;
            if (id != ABPOA_SINK_NODE_ID) fprintf(out_fp, ",");
            else {
                fprintf(out_fp, "\t*\n");
                break;
            }
        }
    }
    free(in_degree);
    for (i = 0; i < n_seq; ++i) free(read_paths[i]); 
    free(read_paths); free(read_path_i);
}

int abpoa_get_aligned_id(abpoa_graph_t *abg, int node_id, uint8_t base) {
    int i, aln_id;
    abpoa_node_t *node = abg->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        aln_id = node[node_id].aligned_node_id[i];
        if (node[aln_id].base == base)
            return aln_id;
    }
    return -1;
}

void abpoa_add_graph_aligned_node1(abpoa_node_t *node, int aligned_id) {
    _uni_realloc(node->aligned_node_id, node->aligned_node_n, node->aligned_node_m, int);
    node->aligned_node_id[node->aligned_node_n++] = aligned_id;
}

void abpoa_add_graph_aligned_node(abpoa_graph_t *abg, int node_id, int aligned_id) {
    int i; abpoa_node_t *node = abg->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        abpoa_add_graph_aligned_node1(node + node[node_id].aligned_node_id[i], aligned_id);
        abpoa_add_graph_aligned_node1(node + aligned_id, node[node_id].aligned_node_id[i]);
    }
    abpoa_add_graph_aligned_node1(abg->node + node_id, aligned_id);
    abpoa_add_graph_aligned_node1(abg->node + aligned_id, node_id);
}

//TODO
void abpoa_set_read_id(uint64_t *read_ids, int read_id) {
    int n = read_id / 64;
    // if ((n = read_id / 64) >= read_ids_n) err_fatal("abpoa_set_read_id", "Unexpected read id: %d\n", read_id);
    uint64_t one = 1; int b = read_id & 0x3f;
    read_ids[n] |= (one << b);
}

int abpoa_add_graph_node(abpoa_graph_t *abg, uint8_t base) {
    int node_id = abg->node_n;
    abpoa_realloc_graph_node(abg);
    // add node
    abg->node[node_id].base = base;
    ++abg->node_n;
    return node_id;
}

int abpoa_add_graph_edge(abpoa_graph_t *abg, int from_id, int to_id, int check_edge, int w, uint8_t add_read_id, int read_id, int read_ids_n) {
    int ret = 1;
    if (from_id < 0 || from_id >= abg->node_n || to_id < 0 || to_id >= abg->node_n) err_fatal(__func__, "node_n: %d\tfrom_id: %d\tto_id: %d.", abg->node_n, from_id, to_id);
    int out_edge_n = abg->node[from_id].out_edge_n;
    if (check_edge) {
        int i;
        for (i = 0; i < out_edge_n; ++i) {
            if (abg->node[from_id].out_id[i] == to_id) { // edge exists
                abg->node[from_id].out_weight[i] += w; // update weight on existing edge
                // update label id
                goto ADD_READ_ID;
                // return;
            }
        }
    }

    // add edge
    /// in edge
    abpoa_realloc_graph_edge(abg, 0, to_id);
    abg->node[to_id].in_id[abg->node[to_id].in_edge_n] = from_id;
    ++abg->node[to_id].in_edge_n;
    /// out edge
    abpoa_realloc_graph_edge(abg, 1, from_id);
    abg->node[from_id].out_id[out_edge_n] = to_id;
    abg->node[from_id].out_weight[out_edge_n] = w; // initial weight for new edge
    ++abg->node[from_id].out_edge_n;
    
    // add read_id to out edge
ADD_READ_ID:
    if (add_read_id) {
        abpoa_node_t *from_node = abg->node + from_id;
        if (from_node->read_ids_n == 0) {
            from_node->read_ids = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t*));
            from_node->read_ids_n = read_ids_n;
        } else if (from_node->read_ids_n < read_ids_n) {
            from_node->read_ids = (uint64_t*)_err_realloc(from_node->read_ids, read_ids_n * sizeof(uint64_t*));
            int i;
            for (i = from_node->read_ids_n; i < read_ids_n; ++i) from_node->read_ids[i] = 0;
            from_node->read_ids_n = read_ids_n;
        }
        abpoa_set_read_id(from_node->read_ids, read_id);
    }
    return ret;
}

void abpoa_add_graph_sequence(abpoa_graph_t *abg, uint8_t *seq, int seq_l, int *qpos_to_node_id, int start, int end, uint8_t add_read_id, int read_id, int read_ids_n) {
    if (start >= seq_l || end <= start) err_fatal(__func__, "seq_l: %d\tstart: %d\tend: %d.", seq_l, start, end);
    if (end > seq_l) end = seq_l;

    int i, last_node_id, cur_node_id;
    last_node_id = ABPOA_SRC_NODE_ID;
    for (i = start; i < end; ++i) {
        cur_node_id = abpoa_add_graph_node(abg, seq[i]);
        if (qpos_to_node_id) qpos_to_node_id[i] = cur_node_id;
        abpoa_add_graph_edge(abg, last_node_id, cur_node_id, 0, 1, add_read_id, read_id, read_ids_n);
        last_node_id = cur_node_id;
    }

    abpoa_add_graph_edge(abg, last_node_id, ABPOA_SINK_NODE_ID, 0, 1, add_read_id, read_id, read_ids_n);
    abg->is_called_cons = abg->is_set_msa_rank = abg->is_topological_sorted = 0;
    // abpoa_topological_sort(abg, abpt);
}

int is_full_upstream_subgraph(abpoa_graph_t *abg, int up_index, int down_index) {
    int i, j, id, in_id;
    for (i = up_index+1; i <= down_index; ++i) {
        id = abg->index_to_node_id[i];
        for (j = 0; j < abg->node[id].in_edge_n; ++j) {
            in_id = abg->node[id].in_id[j];
            if (abg->node_id_to_index[in_id] < up_index) return 0;
        }
    }
    return 1;
}

int abpoa_upstream_index(abpoa_graph_t *abg, int beg_index, int end_index) {
    int min_index, in_index, i, j, node_id, in_id;

    while (1) {
        min_index = beg_index;
        for (i = beg_index; i <= end_index; ++i) {
            node_id = abg->index_to_node_id[i];
            for (j = 0; j < abg->node[node_id].in_edge_n; ++j) {
                in_id = abg->node[node_id].in_id[j];
                in_index = abg->node_id_to_index[in_id];
                min_index = MIN_OF_TWO(min_index, in_index);
            }
        }
        if (is_full_upstream_subgraph(abg, min_index, beg_index)) {
            return min_index;
        } else {
            end_index = beg_index;
            beg_index = min_index; 
        }
    }
}

int is_full_downstream_subgraph(abpoa_graph_t *abg, int up_index, int down_index) {
    int i, j, id, out_id;
    for (i = up_index; i < down_index; ++i) {
        id = abg->index_to_node_id[i];
        for (j = 0; j < abg->node[id].out_edge_n; ++j) {
            out_id = abg->node[id].out_id[j];
            if (abg->node_id_to_index[out_id] > down_index) return 0;
        }
    }
    return 1;
}

int abpoa_downstream_index(abpoa_graph_t *abg, int beg_index, int end_index) {
    int max_index, out_index, i, j, node_id, out_id;

    while (1) {
        max_index = end_index;
        for (i = beg_index; i <= end_index; ++i) {
            node_id = abg->index_to_node_id[i];
            for (j = 0; j < abg->node[node_id].out_edge_n; ++j) {
                out_id = abg->node[node_id].out_id[j];
                out_index = abg->node_id_to_index[out_id];
                max_index = MAX_OF_TWO(max_index, out_index);
            }
        }
        if (is_full_upstream_subgraph(abg, end_index, max_index)) {
            return max_index;
        } else {
            beg_index = end_index;
            end_index = max_index;
        }
    }
}

//   exc_beg | inc_beg ... inc_end | exc_end
void abpoa_subgraph_nodes(abpoa_t *ab, abpoa_para_t *abpt, int exc_beg0, int exc_end0, int *exc_beg, int *exc_end) {
    abpoa_graph_t *abg = ab->abg;
    if (ab->abg->is_topological_sorted == 0) abpoa_topological_sort(abg, abpt);
    int inc_beg_index = abg->node_id_to_index[exc_beg0], inc_end_index = abg->node_id_to_index[exc_end0];

    int exc_beg_index = abpoa_upstream_index(abg, inc_beg_index, inc_end_index);
    int exc_end_index = abpoa_downstream_index(abg, inc_beg_index, inc_end_index);

    if (exc_beg_index < 0 || exc_end_index >= abg->node_n)
        err_fatal_simple("Error in subgraph_nodes");
    *exc_beg = abg->index_to_node_id[exc_beg_index];
    *exc_end = abg->index_to_node_id[exc_end_index];
}

// fusion stratergy :
// 1. Match: merge to one node
// 2. Mismatch: check if B is identical to A' aligned nodes, then merge to node; if not, add node
// 3. Insertion: add node
// 4. Deletion: nothing
// 5. Clipping: add node
// 6. For all first/last node, link to virtual start/end node

// inc_both_ends: set as 1 to add weight for edge between beg_node_id/end_node_id and internal node
int abpoa_add_subgraph_alignment(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *seq, int seq_l, int *qpos_to_node_id, abpoa_res_t res, int read_id, int tot_read_n, int inc_both_ends) {
    abpoa_graph_t *abg = ab->abg;
    int n_cigar = res.n_cigar; abpoa_cigar_t *abpoa_cigar = res.graph_cigar;
    int read_ids_n = 1 + ((tot_read_n-1) >> 6);
    uint8_t add_read_id = abpt->use_read_ids, add, w=1;
    if (abg->node_n == 2) { // empty graph
        abpoa_add_graph_sequence(abg, seq, seq_l, qpos_to_node_id, 0, seq_l, add_read_id, read_id, read_ids_n);
        return 0;
    } else {
        if (abg->node_n < 2) {
            err_fatal(__func__, "Graph node: %d.", abg->node_n);
        } else if (n_cigar == 0) {
            return 0;
            //err_fatal(__func__, "Empty graph cigar.");
        }
    }
    // normal graph, normal graph_cigar
    int i, j; int op, len;
    int node_id, query_id=-1, last_new = 0, last_id = beg_node_id, new_id, aligned_id;

    for (i = 0; i < n_cigar; ++i) {
        op = abpoa_cigar[i] & 0xf;
        if (op == ABPOA_CMATCH) {
            node_id = (abpoa_cigar[i] >> 34) & 0x3fffffff;
            query_id++; // = (abpoa_cigar[i] >> 4) & 0x3fffffff;
            if (abg->node[node_id].base != seq[query_id]) { // mismatch
                // check if query base is identical to node_id's aligned node
                if ((aligned_id = abpoa_get_aligned_id(abg, node_id, seq[query_id])) != -1) {
                    if (last_id != beg_node_id || inc_both_ends) add = 1; else add = 0;
                    abpoa_add_graph_edge(abg, last_id, aligned_id, 1-last_new, w, add_read_id&add, read_id, read_ids_n);
                    last_id = aligned_id; last_new = 0;
                } else {
                    new_id = abpoa_add_graph_node(abg, seq[query_id]);
                    if (last_id != beg_node_id || inc_both_ends) add = 1; else add = 0;
                    abpoa_add_graph_edge(abg, last_id, new_id, 0, w, add_read_id&add, read_id, read_ids_n);
                    last_id = new_id; last_new = 1;
                    // add new_id to node_id's aligned node
                    abpoa_add_graph_aligned_node(abg, node_id, new_id);
                }
            } else { // match
                if (last_id != beg_node_id || inc_both_ends) add = 1; else add = 0;
                abpoa_add_graph_edge(abg, last_id, node_id, 1-last_new, w, add_read_id&add, read_id, read_ids_n);
                last_id = node_id; last_new = 0;
            }
            if (qpos_to_node_id) qpos_to_node_id[query_id] = last_id;
        } else if (op == ABPOA_CINS || op == ABPOA_CSOFT_CLIP || op == ABPOA_CHARD_CLIP) {
            len = (abpoa_cigar[i] >> 4) & 0x3fffffff;
            query_id+=len; // = (abpoa_cigar[i] >> 34) & 0x3fffffff;
            for (j = len-1; j >= 0; --j) { // XXX use dynamic id, instead of static query_id
                new_id = abpoa_add_graph_node(abg, seq[query_id-j]);
                if (last_id != beg_node_id || inc_both_ends) add = 1; else add = 0;
                abpoa_add_graph_edge(abg, last_id, new_id, 0, w, add_read_id&add, read_id, read_ids_n);
                last_id = new_id; last_new = 1;
                if (qpos_to_node_id) qpos_to_node_id[query_id-j] = last_id;
            }
        } else if (op == ABPOA_CDEL) {
            // nothing;
            continue;
        }
    } 
    // if (inc_both_ends) add = 1; else add = 0; XXX end_node_id is always excluded when adding weight
    // abpoa_add_graph_edge(abg, last_id, end_node_id, 1-last_new, w, add_read_id&add, read_id, read_ids_n);
    abpoa_add_graph_edge(abg, last_id, end_node_id, 1-last_new, w, add_read_id, read_id, read_ids_n);
    abg->is_called_cons = abg->is_topological_sorted = 0;
    // abpoa_topological_sort(abg, abpt);
    return 0;
}

int abpoa_add_graph_alignment(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *seq, int seq_l, int *qpos_to_node_id, abpoa_res_t res, int read_id, int tot_read_n, int inc_both_ends) {
    return abpoa_add_subgraph_alignment(ab, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, seq, seq_l, qpos_to_node_id, res, read_id, tot_read_n, inc_both_ends);
}

// reset allocated memery everytime init the graph
// * node
// * index_to_node_id/node_id_to_index/node_id_to_max_remain, max_pos_left/right
void abpoa_reset_graph(abpoa_t *ab, abpoa_para_t *abpt, int qlen) {
    abpoa_graph_t *abg = ab->abg;
    int i, k, node_m;
    abg->is_topological_sorted = abg->is_called_cons = 0;
    for (i = 0; i < abg->node_n; ++i) {
        abg->node[i].in_edge_n = abg->node[i].out_edge_n = abg->node[i].aligned_node_n = 0;
        if (abpt->use_read_ids) {
            for (k = 0; k < abg->node[i].read_ids_n; ++k) 
                abg->node[i].read_ids[k] = 0;
        }
    }
    abg->node_n = 2;
    if (qlen+2 > abg->node_m) {
        node_m = qlen+2; kroundup32(node_m);
        abg->node = (abpoa_node_t*)_err_realloc(abg->node, node_m * sizeof(abpoa_node_t));
        for (i = abg->node_m; i < node_m; ++i) 
            abpoa_set_graph_node(abg, i);
        abg->node_m = abg->index_rank_m = node_m;
        abg->index_to_node_id = (int*)_err_realloc(abg->index_to_node_id, node_m * sizeof(int));
        abg->node_id_to_index = (int*)_err_realloc(abg->node_id_to_index, node_m * sizeof(int));
        if (abpt->out_msa || abpt->cons_agrm == ABPOA_HC || abpt->is_diploid) 
            abg->node_id_to_msa_rank = (int*)_err_realloc(abg->node_id_to_msa_rank, node_m * sizeof(int));
        if (abpt->wb >= 0) {
            abg->node_id_to_max_pos_left = (int*)_err_realloc(abg->node_id_to_max_pos_left, node_m * sizeof(int));
            abg->node_id_to_max_pos_right = (int*)_err_realloc(abg->node_id_to_max_pos_right, node_m * sizeof(int));
            abg->node_id_to_max_remain = (int*)_err_realloc(abg->node_id_to_max_remain, node_m * sizeof(int));
        } else if (abpt->zdrop > 0) {
            abg->node_id_to_max_remain = (int*)_err_realloc(abg->node_id_to_max_remain, node_m * sizeof(int));
        }
    }
    // fprintf(stderr, "qlen: %d, node_n: %d, node_m: %d\n", qlen, abg->node_n, abg->node_m);
}
