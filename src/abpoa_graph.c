#include <stdio.h>
#include <stdlib.h>
#include "abpoa_align.h"
#include "simd_abpoa_align.h"
#include "kdq.h"

extern char LogTable65536[65536];
extern char bit_table16[65536];

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

static inline int ilog2_64(uint64_t v)
{
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

void abpoa_set_graph_node(abpoa_graph_t *graph, int node_i) {
    graph->node[node_i].node_id = node_i;
    graph->node[node_i].in_edge_n = 0; graph->node[node_i].in_edge_m = 0;
    graph->node[node_i].out_edge_n = 0; graph->node[node_i].out_edge_m = 0;
    graph->node[node_i].aligned_node_n = 0; graph->node[node_i].aligned_node_m = 0;
    graph->node[node_i].read_ids_n = 0;
}

void abpoa_free_node(abpoa_para_t *abpt, abpoa_node_t *node, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (node[i].in_edge_m > 0) free(node[i].in_id);
        if (node[i].out_edge_m > 0) { 
            free(node[i].out_id); free(node[i].out_weight);
            if (abpt->use_read_ids) {
                if (node[i].read_ids_n > 0)
                    free(node[i].read_ids);
            }
        }
        if (node[i].aligned_node_m > 0) free(node[i].aligned_node_id);
    }
    free(node);
}

// 0: in_edge, 1: out_edge
abpoa_graph_t *abpoa_realloc_graph_edge(abpoa_graph_t *graph, int io, int id) {
    if (io == 0) {
        _uni_realloc(graph->node[id].in_id, graph->node[id].in_edge_n, graph->node[id].in_edge_m, int);
    } else {
        int edge_m = graph->node[id].out_edge_m;
        _uni_realloc(graph->node[id].out_id, graph->node[id].out_edge_n, edge_m, int);
        _uni_realloc(graph->node[id].out_weight, graph->node[id].out_edge_n, graph->node[id].out_edge_m, int);
    }
    return graph;
}

abpoa_graph_t *abpoa_realloc_graph_node(abpoa_graph_t *graph) {
    if (graph->node_m <= 0) {
        graph->node_m = 1;
        graph->node = (abpoa_node_t*)_err_calloc(1, sizeof(abpoa_node_t));
    }
    if (graph->node_n == graph->node_m) {
        int i;
        graph->node_m <<= 1;
        graph->node = (abpoa_node_t*)_err_realloc(graph->node, graph->node_m * sizeof(abpoa_node_t));
        for (i = graph->node_m >> 1; i < graph->node_m; ++i) {
            abpoa_set_graph_node(graph, i);
        }
    }
    return graph;
}

abpoa_graph_t *abpoa_init_graph(void) {
    abpoa_graph_t *graph = (abpoa_graph_t*)_err_malloc(sizeof(abpoa_graph_t));
    graph->node_n = 2, graph->node_m = 2, graph->index_rank_m = 0;
    graph->node = abpoa_init_node(2);
    graph->node[0].node_id = 0; graph->node[1].node_id = 1;
    graph->node[0].read_ids_n = 0; graph->node[1].read_ids_n = 0;
    graph->cons_l = graph->cons_m = 0; graph->cons_seq = NULL;
    graph->is_topological_sorted = graph->is_called_cons = 0;
    graph->node_id_to_index = NULL; graph->index_to_node_id = NULL; graph->node_id_to_msa_rank = NULL;
    graph->node_id_to_max_pos_left = NULL; graph->node_id_to_max_pos_right = NULL; graph->node_id_to_max_remain = NULL;
    return graph;
}

void abpoa_free_graph(abpoa_graph_t *graph, abpoa_para_t *abpt) {
    if (graph->node_m > 0) abpoa_free_node(abpt, graph->node, graph->node_m);
    if (graph->cons_m > 0) free(graph->cons_seq);

    if (graph->node_n > 0) {
        free(graph->index_to_node_id);
        free(graph->node_id_to_index);
        if (graph->node_id_to_msa_rank) free(graph->node_id_to_msa_rank);

        if (graph->node_id_to_max_pos_left) free(graph->node_id_to_max_pos_left);
        if (graph->node_id_to_max_pos_right) free(graph->node_id_to_max_pos_right);
        if (graph->node_id_to_max_remain) free(graph->node_id_to_max_remain);
    }
    free(graph);
}

abpoa_t *abpoa_init(void) {
    abpoa_t *ab = (abpoa_t*)_err_malloc(sizeof(abpoa_t));
    ab->abg = abpoa_init_graph();
    ab->abm = abpoa_init_simd_matrix();
    return ab;
}

void abpoa_free(abpoa_t *ab, abpoa_para_t *abpt) {
    abpoa_free_graph(ab->abg, abpt);
    abpoa_free_simd_matrix(ab->abm);
    free(ab);
}

void abpoa_BFS_set_node_index(abpoa_graph_t *graph, int src_id, int sink_id) {
    int *id, cur_id, i, j, out_id, aligned_id, q_size, new_q_size;
    int index = 0;

    int *in_degree = (int*)_err_malloc(graph->node_n * sizeof(int));
    for (i = 0; i < graph->node_n; ++i) in_degree[i] = graph->node[i].in_edge_n;

    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, src_id); q_size = 1; new_q_size = 0; // node[q.id].in_degree equals 0
    while (q_size > 0) {
        if ((id = kdq_shift_int(q)) == 0) err_fatal_simple("Error in queue.");
        cur_id = *id;
        graph->index_to_node_id[index] = cur_id;
        graph->node_id_to_index[cur_id] = index++;

        if (cur_id == sink_id) {
            kdq_destroy_int(q); free(in_degree);
            return;
        }
        for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
            out_id = graph->node[cur_id].out_id[i];
            if (--in_degree[out_id] == 0) {
                for (j = 0; j < graph->node[out_id].aligned_node_n; ++j) {
                    aligned_id = graph->node[out_id].aligned_node_id[j];
                    if (in_degree[aligned_id] != 0) goto next_out_node;
                }
                kdq_push_int(q, out_id);
                ++new_q_size;
                for (j = 0; j < graph->node[out_id].aligned_node_n; ++j) {
                    aligned_id = graph->node[out_id].aligned_node_id[j];
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

void abpoa_BFS_set_node_remain(abpoa_graph_t *graph, int src_id, int sink_id) {
    int *id, cur_id, i, out_id, in_id;

    int *out_degree = (int*)_err_malloc(graph->node_n * sizeof(int));
    for (i = 0; i < graph->node_n; ++i) {
        out_degree[i] = graph->node[i].out_edge_n;
        graph->node_id_to_max_remain[i] = 0;
    }

    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, sink_id); // node[q.id].in_degree equals 0
    graph->node_id_to_max_remain[sink_id] = -1; // XXX not 0
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;

        // all out_id of cur_id have beed visited
        // max weight out_id
        if (cur_id != sink_id) {
            int max_w=-1, max_id=-1;
            for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
                out_id = graph->node[cur_id].out_id[i];
                if (graph->node[cur_id].out_weight[i] > max_w) {
                    max_w = graph->node[cur_id].out_weight[i];
                    max_id = out_id;
                }
            }
            graph->node_id_to_max_remain[cur_id] = graph->node_id_to_max_remain[max_id] + 1;
            // fprintf(stderr, "%d -> %d\n", graph->node_id_to_index[cur_id], graph->node_id_to_max_remain[cur_id]);
        }
        if (cur_id == src_id) {
            kdq_destroy_int(q); free(out_degree);
            return;
        }
        for (i = 0; i < graph->node[cur_id].in_edge_n; ++i) {
            in_id = graph->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) kdq_push_int(q, in_id);
        }
    }
    err_fatal_simple("Failed to set node remain.");
}

// 1. index_to_node_id
// 2. node_id_to_index
// 3. node_id_to_rank
void abpoa_topological_sort(abpoa_t *ab, abpoa_para_t *abpt) {
    abpoa_graph_t *graph = ab->abg;
    if (graph->node_n <= 0) {
        err_func_format_printf(__func__, "Empty graph.\n");
        return;
    }
    int node_n = graph->node_n;
    if (node_n > graph->index_rank_m) {
        graph->index_rank_m = node_n; kroundup32(graph->index_rank_m);
        // fprintf(stderr, "node_n: %d, index_rank_m: %d\n", node_n, graph->index_rank_m);
        graph->index_to_node_id = (int*)_err_realloc(graph->index_to_node_id, graph->index_rank_m * sizeof(int));
        graph->node_id_to_index = (int*)_err_realloc(graph->node_id_to_index, graph->index_rank_m * sizeof(int));
        if (abpt->out_msa || abpt->cons_agrm == ABPOA_HC || abpt->is_diploid) 
            graph->node_id_to_msa_rank = (int*)_err_realloc(graph->node_id_to_msa_rank, graph->index_rank_m * sizeof(int));
        if (abpt->wb >= 0) {
            graph->node_id_to_max_pos_left = (int*)_err_realloc(graph->node_id_to_max_pos_left, graph->index_rank_m * sizeof(int));
            graph->node_id_to_max_pos_right = (int*)_err_realloc(graph->node_id_to_max_pos_right, graph->index_rank_m * sizeof(int));
            graph->node_id_to_max_remain = (int*)_err_realloc(graph->node_id_to_max_remain, graph->index_rank_m * sizeof(int));
        } else if (abpt->zdrop > 0) {
            graph->node_id_to_max_remain = (int*)_err_realloc(graph->node_id_to_max_remain, graph->index_rank_m * sizeof(int));
        }
    }
    // start from ABPOA_SRC_NODE_ID to ABPOA_SINK_NODE_ID
    abpoa_BFS_set_node_index(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    // init min/max rank
    if (abpt->wb >= 0) {
        int i;
        for (i = 0; i < node_n; ++i) {
            graph->node_id_to_max_pos_right[i] = 0;
            graph->node_id_to_max_pos_left[i] = node_n;
        }
        abpoa_BFS_set_node_remain(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    } else if (abpt->zdrop > 0)
        abpoa_BFS_set_node_remain(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    graph->is_topological_sorted = 1;
}

void abpoa_DFS_set_msa_rank(abpoa_graph_t *graph, int src_id, int sink_id, int *in_degree) {
    // fprintf(stderr, "node_n: %d, m: %d\n", graph->node_n, graph->index_rank_m);
    if (graph->node_n > graph->index_rank_m) {
        int m = graph->node_n; kroundup32(m);
        graph->node_id_to_msa_rank = (int*)_err_realloc(graph->node_id_to_msa_rank, m * sizeof(int));
    }
    int *id, cur_id, i, j, out_id, aligned_id;
    int msa_rank = 0;
    kdq_int_t *q = kdq_init_int();

    // Depth-First-Search
    kdq_push_int(q, src_id); // node[q.id].in_degree equals 0
    graph->node_id_to_msa_rank[src_id] = -1;
    // printf("tot_node_n: %d, node_m: %d\n", graph->node_n, graph->node_m);

    while((id = kdq_pop_int(q)) != 0) {
        cur_id = *id;
        if (graph->node_id_to_msa_rank[cur_id] < 0) {
            graph->node_id_to_msa_rank[cur_id] = msa_rank;
            for (i = 0; i < graph->node[cur_id].aligned_node_n; ++i) {
                aligned_id = graph->node[cur_id].aligned_node_id[i];
                graph->node_id_to_msa_rank[aligned_id] = msa_rank;
            }
            msa_rank++;
        }

        if (cur_id == sink_id) {
            kdq_destroy_int(q);
            graph->is_set_msa_rank = 1;
            return;
        }
        for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
            out_id = graph->node[cur_id].out_id[i];
            if (--in_degree[out_id] == 0) {
                for (j = 0; j < graph->node[out_id].aligned_node_n; ++j) {
                    aligned_id = graph->node[out_id].aligned_node_id[j];
                    if (in_degree[aligned_id] != 0) goto next_out_node;
                }
                kdq_push_int(q, out_id);
                graph->node_id_to_msa_rank[out_id] = -1;
                for (j = 0; j < graph->node[out_id].aligned_node_n; ++j) {
                    aligned_id = graph->node[out_id].aligned_node_id[j];
                    kdq_push_int(q, aligned_id);
                    // printf("aln_id: %d\n", aligned_id);
                    graph->node_id_to_msa_rank[aligned_id] = -1;
                }
            }
next_out_node:;
        }
    }
    err_fatal_simple("Error in set_msa_rank.\n");
}

void abpoa_set_msa_rank(abpoa_graph_t *graph, int src_id, int sink_id) {
    if (graph->is_set_msa_rank == 0) {
        int i, *in_degree = (int*)_err_malloc(graph->node_n * sizeof(int));
        for (i = 0; i < graph->node_n; ++i) in_degree[i] = graph->node[i].in_edge_n;
        abpoa_DFS_set_msa_rank(graph, src_id, sink_id, in_degree);
        free(in_degree);
    }
}

void abpoa_set_row_column_weight(abpoa_graph_t *graph, int **rc_weight) {
    int i, k, rank, aligned_id;
    uint64_t b;
    for (i = 2; i < graph->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(graph, i);
        for (k = 0; k < graph->node[i].aligned_node_n; ++k) {
            aligned_id = graph->node[i].aligned_node_id[k];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(graph, aligned_id));
        }
        // assign seq
        for (k = 0; k < graph->node[i].read_ids_n; ++k) {
            b = graph->node[i].read_ids[k];
            rc_weight[rank-1][graph->node[i].base] += get_bit_cnt4(bit_table16, b);
        }
        rc_weight[rank-1][4] -= rc_weight[rank-1][graph->node[i].base];
    }
}

void abpoa_set_row_column_ids_weight(abpoa_graph_t *graph, uint64_t ***read_ids, int **rc_weight, int msa_l, int seq_n, int read_ids_n) {
    int i, j, k, rank, aligned_id;
    uint64_t b, one = 1, *whole_read_ids = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    for (i = 0; i < seq_n; ++i) {
        j = i / 64; b = i & 0x3f;
        whole_read_ids[j] |= (one << b);
    }
    for (i = 0; i < msa_l; ++i) {
        for (j = 0; j < read_ids_n; ++j) {
            read_ids[i][4][j] = whole_read_ids[j];
        }
    }
    for (i = 2; i < graph->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(graph, i);
        for (k = 0; k < graph->node[i].aligned_node_n; ++k) {
            aligned_id = graph->node[i].aligned_node_id[k];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(graph, aligned_id));
        }
        // assign seq
        for (k = 0; k < graph->node[i].read_ids_n; ++k) {
            b = graph->node[i].read_ids[k];
            rc_weight[rank-1][graph->node[i].base] += get_bit_cnt4(bit_table16, b);
            read_ids[rank-1][graph->node[i].base][k] = b;
            read_ids[rank-1][4][k] ^= b;
        }
        rc_weight[rank-1][4] -= rc_weight[rank-1][graph->node[i].base];
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

void output_consensus(abpoa_graph_t *graph, FILE *out_fp) {
    int i;
    fprintf(out_fp, ">Consensus_sequence\n");
    for (i = 0; i < graph->cons_l; ++i) {
        fprintf(out_fp, "%c", "ACGTN"[graph->cons_seq[i]]);
    } fprintf(out_fp, "\n");
}

void abpoa_store_consensus(abpoa_graph_t *graph, uint8_t ***cons_seq, int **cons_l) {
    *cons_seq = (uint8_t**)_err_malloc(sizeof(uint8_t*));
    *cons_l = (int*)_err_malloc(sizeof(int));
    (*cons_seq)[0] = (uint8_t*)_err_malloc(sizeof(uint8_t) * graph->cons_l);
    (*cons_l)[0] = graph->cons_l;
    int i;
    for (i = 0; i < graph->cons_l; ++i) (*cons_seq)[0][i] = graph->cons_seq[i];
}

void abpoa_generate_consensus_core(abpoa_graph_t *graph, int path_start, int sink_id, int *max_out_id) {
    if (graph->node_n > graph->cons_m) {
        graph->cons_m = MAX_OF_TWO(graph->cons_m << 1, graph->node_n);
        graph->cons_seq = (uint8_t*)_err_realloc(graph->cons_seq, graph->cons_m * sizeof(uint8_t));
    }
    int id = path_start;
    while (id != sink_id) {
        graph->cons_seq[graph->cons_l++] = graph->node[id].base;
        id = max_out_id[id];
    }
}

void abpoa_traverse_min_flow(abpoa_graph_t *graph,  int src_id, int sink_id, int *out_degree) {
    int *id, i, cur_id, in_id, out_id, max_id, max_w, max_out_i, min_w;
    kdq_int_t *q = kdq_init_int(); kdq_push_int(q, sink_id);
    int *max_out_id = (int*)_err_malloc(graph->node_n * sizeof(int)), *max_out_weight = (int*)_err_malloc(graph->node_n * sizeof(int));
    int path_start=-1;
    // reverse Breadth-First-Search
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == sink_id) {
            max_out_id[sink_id] = -1;
            max_out_weight[sink_id] = INT32_MAX;
        } else {
            max_w = INT32_MIN, max_id = -1, max_out_i = -1, min_w = INT32_MAX;
            for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
                out_id = graph->node[cur_id].out_id[i];
                min_w = MIN_OF_TWO(max_out_weight[out_id], graph->node[cur_id].out_weight[i]);
                if (max_w < min_w) {
                    max_w = min_w;
                    max_id = out_id;
                    max_out_i = i;
                } else if (max_w == min_w && graph->node[cur_id].out_weight[max_out_i] <= graph->node[cur_id].out_weight[i]) {
                    max_id = out_id;
                    max_out_i = i;
                }
            }
            max_out_id[cur_id] = max_id;
            max_out_weight[cur_id] = max_w;
        }

        if (cur_id == src_id) { 
            path_start = max_out_id[cur_id];
            kdq_destroy_int(q); 
            goto MF_CONS; 
        }
        for (i = 0; i < graph->node[cur_id].in_edge_n; ++i) {
            in_id = graph->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) kdq_push_int(q, in_id);
        }
    }
MF_CONS:
    abpoa_generate_consensus_core(graph, path_start, sink_id, max_out_id);
    free(max_out_id); free(max_out_weight);
}

// heaviest_bundling
// 1. argmax{cur->weight}
// 2. argmax{out_node->weight}
void abpoa_heaviest_bundling(abpoa_graph_t *graph,  int src_id, int sink_id, int *out_degree) {
    int *id, i, cur_id, in_id, out_id, max_id, max_w, out_w;
    int *max_out_id = (int*)_err_malloc(graph->node_n * sizeof(int)), *score = (int*)_err_malloc(graph->node_n * sizeof(int));
    int path_start = -1, path_score = -1;
    kdq_int_t *q = kdq_init_int();
    kdq_push_int(q, sink_id);
    // reverse Breadth-First-Search
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == sink_id) {
            max_out_id[sink_id] = -1;
            score[sink_id] = 0;
        } else {
            max_w = INT32_MIN, max_id = -1;
            for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
                out_id = graph->node[cur_id].out_id[i];
                out_w = graph->node[cur_id].out_weight[i];
                if (max_w < out_w) {
                    max_w = out_w;
                    max_id = out_id;
                } else if (max_w == out_w && score[max_id] <= score[out_id]) {
                    max_id = out_id;
                }
            }
            score[cur_id] = max_w + score[max_id];
            max_out_id[cur_id] = max_id;
        }
        if (cur_id == src_id) {
            for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
                out_id = graph->node[cur_id].out_id[i];
                if (score[out_id] > path_score) {
                    path_start = out_id;
                    path_score = score[out_id];
                }
            }
            kdq_destroy_int(q);
            goto HB_CONS;
        }
        for (i = 0; i < graph->node[cur_id].in_edge_n; ++i) {
            in_id = graph->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) 
                kdq_push_int(q, in_id);
        }
    }
HB_CONS:
    abpoa_generate_consensus_core(graph, path_start, sink_id, max_out_id);
    free(max_out_id); free(score);
}

void abpoa_heaviest_column_consensus(abpoa_graph_t *graph, int **rc_weight, int msa_l, int seq_n) {
    if (graph->cons_m < msa_l+1) {
        graph->cons_m = msa_l; kroundup32(graph->cons_m);
        graph->cons_seq = (uint8_t*)_err_realloc(graph->cons_seq, graph->cons_m * sizeof(uint8_t));
    }
    int i, j, w, max_base, max_w, gap_w;
    for (i = 0; i < msa_l; ++i) {
        max_w = 0, max_base = 5, gap_w = seq_n;
        for (j = 0; j < 4; ++j) {
            w = rc_weight[i][j];
            if (w > max_w) {
                max_base = j, max_w = w;
            }
            gap_w -= w;
        }
        if (max_w > gap_w) {
            graph->cons_seq[graph->cons_l++] = max_base;
        }
    }
}

void abpoa_heaviest_column(abpoa_graph_t *graph, int src_id, int sink_id, int seq_n) {
    abpoa_set_msa_rank(graph, src_id, sink_id);

    int i, msa_l = graph->node_id_to_msa_rank[sink_id] - 1;
    int **rc_weight = (int**)_err_malloc(sizeof(int*) * msa_l);
    for (i = 0; i < msa_l; ++i) {
        rc_weight[i] = (int*)_err_calloc(5, sizeof(int)); // ACGT
        rc_weight[i][4] = seq_n;
    } 
    abpoa_set_row_column_weight(graph, rc_weight);

    abpoa_heaviest_column_consensus(graph, rc_weight, msa_l, seq_n);
    for (i = 0; i < msa_l; ++i) free(rc_weight[i]); free(rc_weight);
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

int set_clu_read_ids(int **clu_read_ids, int *clu_read_ids_n, int seq_n, double min_freq, uint64_t ***read_ids, uint8_t **het_cons_base, int *het_pos, int het_n) {
    int i, j, seq_i, pos, clu_base;
    int **diff = (int**)_err_malloc(sizeof(int*) * seq_n);
    for (i = 0; i < seq_n; ++i) {
        diff[i] = (int*)_err_malloc(2 *  sizeof(int));
        diff[i][0] = het_n, diff[i][1] = het_n;
    }
    uint64_t seq_b, one = 1; int seq_bit;
    for (seq_i = 0; seq_i < seq_n; ++seq_i) {
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
    for (i = 0; i < seq_n; ++i) {
        if (diff[i][0] < diff[i][1]) {
            clu_read_ids[0][clu_read_ids_n[0]++] = i;
        } else if (diff[i][1] < diff[i][0]) {
            clu_read_ids[1][clu_read_ids_n[1]++] = i;
        }
    }
    for (i = 0; i < seq_n; ++i) free(diff[i]); free(diff);

    int min_n = MIN_OF_TWO(clu_read_ids_n[0], clu_read_ids_n[1]);
    int clu_n = (min_n >= (int)(seq_n * min_freq)) ? 2 : 1; 
    return clu_n;
}

int abpoa_diploid_ids(uint64_t ***read_ids, int **rc_weight, int msa_l, int seq_n, double min_freq, int read_ids_n, int **clu_read_ids, int *clu_read_ids_n) {
    int i, j;
    uint64_t **het_read_ids = (uint64_t**)_err_malloc(sizeof(uint64_t*) * 2 );
    for (i = 0; i < 2; ++i) het_read_ids[i] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    int *het_pos = (int*)_err_malloc(sizeof(int) * msa_l), het_n = 0;
    uint8_t **het_cons_base = (uint8_t**)_err_malloc(sizeof(uint8_t*) * msa_l);
    for (i = 0; i < msa_l; ++i) het_cons_base[i] = (uint8_t*)_err_malloc(sizeof(uint8_t) * 2);
    int min_w = MAX_OF_TWO(1, seq_n * min_freq);
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
        clu_n = set_clu_read_ids(clu_read_ids, clu_read_ids_n, seq_n, min_freq, read_ids, het_cons_base, het_pos, het_n);
    }
    
    for (i = 0; i < msa_l; ++i) free(het_cons_base[i]);
    for (i = 0; i < 2; ++i) free(het_read_ids[i]); 
    free(het_cons_base); free(het_read_ids); free(het_pos);
    return clu_n;
}

void abpoa_diploid_heaviest_column(abpoa_graph_t *graph, int src_id, int sink_id, int seq_n, double min_freq, FILE *out_fp, uint8_t ***cons_seq, int **cons_l, int *cons_n) {
    abpoa_set_msa_rank(graph, src_id, sink_id);
    int i, j, msa_l = graph->node_id_to_msa_rank[sink_id] - 1;
    int read_ids_n = (seq_n-1)/64+1;
    uint64_t ***read_ids = (uint64_t***)_err_malloc(sizeof(uint64_t**) * msa_l);
    for (i = 0; i < msa_l; ++i) {
        read_ids[i] = (uint64_t**)_err_malloc(sizeof(uint64_t*) * 5);
        for (j = 0; j < 5; ++j) read_ids[i][j] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    }
    int **rc_weight = (int**)_err_malloc(msa_l * sizeof(int*));
    for (i = 0; i < msa_l; ++i) {
        rc_weight[i] = (int*)_err_calloc(5, sizeof(int)); // ACGT
        rc_weight[i][4] = seq_n;
    } 
    abpoa_set_row_column_ids_weight(graph, read_ids, rc_weight, msa_l, seq_n, read_ids_n);

    int **clu_read_ids = (int**)_err_malloc(sizeof(int*) * 2), *clu_read_ids_n;
    for (i = 0; i < 2; ++i) clu_read_ids[i] = (int*)_err_malloc(sizeof(int) * seq_n);
    clu_read_ids_n = (int*)_err_calloc(2, sizeof(int));
    int clu_n = abpoa_diploid_ids(read_ids, rc_weight, msa_l, seq_n, min_freq, read_ids_n, clu_read_ids, clu_read_ids_n);
    if (clu_n == 1) {
        abpoa_heaviest_column_consensus(graph, rc_weight, msa_l, seq_n);
        if (out_fp) output_consensus(graph, out_fp);
        if (cons_n) {
            *cons_n = 1; abpoa_store_consensus(graph, cons_seq, cons_l);
        }
    } else abpoa_heaviest_column_multip_consensus(read_ids, clu_read_ids, clu_read_ids_n, clu_n, read_ids_n, msa_l, out_fp, cons_seq, cons_l, cons_n);

    for (i = 0; i < msa_l; ++i) {
        for (j = 0; j < 5; ++j) {
            free(read_ids[i][j]);
        } free(read_ids[i]); free(rc_weight[i]);
    } free(read_ids); free(rc_weight);
    for (i = 0; i < 2; ++i) free(clu_read_ids[i]); free(clu_read_ids); free(clu_read_ids_n);
}

// should always topological sort first, then generate consensus
int abpoa_generate_consensus(abpoa_t *ab, abpoa_para_t *abpt, int seq_n, FILE *out_fp, uint8_t ***cons_seq, int **cons_l, int *cons_n) {
    abpoa_graph_t *graph = ab->abg;
    if (graph->node_n <= 2) return 0;
    int i, *out_degree = (int*)_err_malloc(graph->node_n * sizeof(int));
    for (i = 0; i < graph->node_n; ++i) {
        out_degree[i] = graph->node[i].out_edge_n;
    }

    if (abpt->is_diploid) {
        abpoa_diploid_heaviest_column(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, seq_n, abpt->min_freq, out_fp, cons_seq, cons_l, cons_n);
    } else {
        if (abpt->cons_agrm == ABPOA_HB) abpoa_heaviest_bundling(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree);
        else if (abpt->cons_agrm == ABPOA_HC) abpoa_heaviest_column(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, seq_n);
        else if (abpt->cons_agrm == ABPOA_MF) abpoa_traverse_min_flow(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree); 
        else err_fatal(__func__, "Unknown consensus calling algorithm: %d.", abpt->cons_agrm);
        if (out_fp) output_consensus(graph, out_fp);
        if (cons_n) {
            *cons_n = 1; abpoa_store_consensus(graph, cons_seq, cons_l);
        }
    } 
    graph->is_called_cons = 1; free(out_degree);
    return graph->cons_l;
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
            // printf("%d -> %d\n", node.node_id, read_id);
            msa_seq[b+read_id][rank-1] = base;
            num ^= tmp;
        }
        b += 64;
    }
}

void abpoa_generate_rc_msa(abpoa_t *ab, int seq_n, FILE *out_fp, uint8_t ***msa_seq, int *msa_l) {
    abpoa_graph_t *graph = ab->abg;
    if (graph->node_n <= 2) return;
    abpoa_set_msa_rank(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);

    int i, j, k, aligned_id, rank;
    uint8_t **_msa_seq = (uint8_t**)_err_malloc(seq_n * sizeof(uint8_t*));
    int _msa_l = graph->node_id_to_msa_rank[ABPOA_SINK_NODE_ID]-1;

    for (i = 0; i < seq_n; ++i) {
        _msa_seq[i] = (uint8_t*)_err_malloc(_msa_l * sizeof(uint8_t));
        for (j = 0; j < _msa_l; ++j) 
            _msa_seq[i][j] = 5;
    }

    if (out_fp) fprintf(out_fp, ">Multiple_sequence_alignment\n");
    for (i = 2; i < graph->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(graph, i);
        for (k = 0; k < graph->node[i].aligned_node_n; ++k) {
            aligned_id = graph->node[i].aligned_node_id[k];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(graph, aligned_id));
        }
        // assign seq
        abpoa_set_msa_seq(graph->node[i], rank, _msa_seq);
    }
    if (out_fp) {
        for (i = 0; i < seq_n; ++i) {
            for (j = 0; j < _msa_l; ++j) fprintf(out_fp, "%c", "ACGTN-"[_msa_seq[i][j]]);
            fprintf(out_fp, "\n");
        }
    }
    if (msa_l) {
        *msa_l = _msa_l; *msa_seq = _msa_seq;
    } else {
        for (i = 0; i < seq_n; ++i) free(_msa_seq[i]); free(_msa_seq);
    }
}

int abpoa_get_aligned_id(abpoa_graph_t *graph, int node_id, uint8_t base) {
    int i;
    abpoa_node_t *node = graph->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        if (node[node[node_id].aligned_node_id[i]].base == base)
            return graph->node[node_id].aligned_node_id[i];
    }
    return -1;
}

void abpoa_add_graph_aligned_node1(abpoa_node_t *node, int aligned_id) {
    _uni_realloc(node->aligned_node_id, node->aligned_node_n, node->aligned_node_m, int);
    node->aligned_node_id[node->aligned_node_n++] = aligned_id;
}

void abpoa_add_graph_aligned_node(abpoa_t *ab, int node_id, int aligned_id) {
    abpoa_graph_t *graph = ab->abg;
    int i; abpoa_node_t *node = graph->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        abpoa_add_graph_aligned_node1(node + node[node_id].aligned_node_id[i], aligned_id);
        abpoa_add_graph_aligned_node1(node + aligned_id, node[node_id].aligned_node_id[i]);
    }
    abpoa_add_graph_aligned_node1(graph->node + node_id, aligned_id);
    abpoa_add_graph_aligned_node1(graph->node + aligned_id, node_id);
}

//TODO
void abpoa_set_read_id(uint64_t *read_ids, int read_id) {
    int n = read_id / 64;
    // if ((n = read_id / 64) >= read_ids_n) err_fatal("abpoa_set_read_id", "Unexpected read id: %d\n", read_id);
    uint64_t one = 1; int b = read_id & 0x3f;
    read_ids[n] |= (one << b);
}

int abpoa_add_graph_node(abpoa_t *ab, uint8_t base) {
    abpoa_graph_t *graph = ab->abg;
    int node_id = graph->node_n;
    graph = abpoa_realloc_graph_node(graph);
    // add node
    graph->node[node_id].base = base;
    ++graph->node_n;
    return node_id;
}

int abpoa_add_graph_edge(abpoa_t *ab, int from_id, int to_id, int check_edge, uint8_t add_read_id, int read_id, int read_ids_n) {
    int ret = 1;
    abpoa_graph_t *graph = ab->abg;
    if (from_id < 0 || from_id >= graph->node_n || to_id < 0 || to_id >= graph->node_n) err_fatal(__func__, "node_n: %d\tfrom_id: %d\tto_id: %d.", graph->node_n, from_id, to_id);
    int out_edge_n = graph->node[from_id].out_edge_n;
    if (check_edge) {
        int i;
        for (i = 0; i < out_edge_n; ++i) {
            if (graph->node[from_id].out_id[i] == to_id) { // edge exists
                graph->node[from_id].out_weight[i]++; // update weight on existing edge
                // update label id
                goto ADD_READ_ID;
                // return;
            }
        }
    }

    // add edge
    /// in edge
    graph = abpoa_realloc_graph_edge(graph, 0, to_id);
    graph->node[to_id].in_id[graph->node[to_id].in_edge_n] = from_id;
    ++graph->node[to_id].in_edge_n;
    /// out edge
    graph = abpoa_realloc_graph_edge(graph, 1, from_id);
    graph->node[from_id].out_id[out_edge_n] = to_id;
    graph->node[from_id].out_weight[out_edge_n] = 1; // initial weight for new edge
    ++graph->node[from_id].out_edge_n;
    
    // add read_id to out edge
ADD_READ_ID:
    if (add_read_id) {
        abpoa_node_t *from_node = graph->node + from_id;
        if (from_node->read_ids_n == 0) {
            from_node->read_ids = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t*));
            from_node->read_ids_n = read_ids_n;
        } else if (from_node->read_ids_n < read_ids_n) {
            from_node->read_ids = (uint64_t*)_err_realloc(from_node->read_ids, read_ids_n * sizeof(uint64_t*));
            from_node->read_ids_n = read_ids_n;
        }
        abpoa_set_read_id(from_node->read_ids, read_id);
    }
    return ret;
}

void abpoa_add_graph_sequence(abpoa_t *ab, uint8_t *seq, int seq_l, int start, int end, uint8_t add_read_id, int read_id, int read_ids_n) {
    if (seq_l <= 0 || start >= seq_l || end <= start) err_fatal(__func__, "seq_l: %d\tstart: %d\tend: %d.", seq_l, start, end);
    if (start < 0) start = 0; if (end > seq_l) end = seq_l;
    abpoa_graph_t *graph = ab->abg;
    int node_id = abpoa_add_graph_node(ab, seq[start]);
    abpoa_add_graph_edge(ab, ABPOA_SRC_NODE_ID, node_id, 0, add_read_id, read_id, read_ids_n);
    int i; 
    for (i = start+1; i < end; ++i) {
        node_id = abpoa_add_graph_node(ab, seq[i]);
        abpoa_add_graph_edge(ab, node_id-1, node_id, 0, add_read_id, read_id, read_ids_n);
    }
    abpoa_add_graph_edge(ab, node_id, ABPOA_SINK_NODE_ID, 0, add_read_id, read_id, read_ids_n);
    graph->is_called_cons = graph->is_set_msa_rank = graph->is_topological_sorted = 0;
    // abpoa_topological_sort(graph, abpt);
}

// fusion stratergy :
// 1. Match: merge to one node
// 2. Mismatch: check if B is identical to A' aligned nodes, then merge to node; if not, add node
// 3. Insertion: add node
// 4. Deletion: nothing
// 5. Clipping: add node
// 6. For all first/last node, link to virtual start/end node
//33S:32
//1M:74,33
//26I:59

// for seq-to-graph alignment
//     generate dot file
int abpoa_add_graph_alignment(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *seq, int seq_l, int n_cigar, abpoa_cigar_t *abpoa_cigar, int read_id, int tot_read_n) {
    abpoa_graph_t *graph = ab->abg; int read_ids_n = 1 + ((tot_read_n-1) >> 6);
    uint8_t add_read_id = abpt->use_read_ids;
    if (graph->node_n == 2) { // empty graph
        abpoa_add_graph_sequence(ab, seq, seq_l, 0, seq_l, add_read_id, read_id, read_ids_n);
        return 0;
    } else {
        if (graph->node_n < 2) {
            err_fatal(__func__, "Graph node: %d.", graph->node_n);
        } else if (n_cigar == 0) {
            return 0;
            //err_fatal(__func__, "Empty graph cigar.");
        }
    }
    // normal graph, normal graph_cigar
    int i, j; int op, len, node_id, query_id, last_new = 0, last_id = ABPOA_SRC_NODE_ID, new_id, aligned_id;
    for (i = 0; i < n_cigar; ++i) {
        op = abpoa_cigar[i] & 0xf;
        if (op == ABPOA_CMATCH) {
            node_id = (abpoa_cigar[i] >> 34) & 0x3fffffff;
            query_id = (abpoa_cigar[i] >> 4) & 0x3fffffff;
            if (graph->node[node_id].base != seq[query_id]) { // mismatch
                // check if query base is identical to node_id's aligned node
                if ((aligned_id = abpoa_get_aligned_id(graph, node_id, seq[query_id])) >= 0) {
                    abpoa_add_graph_edge(ab, last_id, aligned_id, 1-last_new, add_read_id, read_id, read_ids_n);
                    last_id = aligned_id; last_new = 0;
                } else {
                    new_id = abpoa_add_graph_node(ab, seq[query_id]);
                    abpoa_add_graph_edge(ab, last_id, new_id, 0, add_read_id, read_id, read_ids_n);
                    last_id = new_id; last_new = 1;
                    // add new_id to node_id's aligned node
                    abpoa_add_graph_aligned_node(ab, node_id, new_id);
                }
            } else { // match
                abpoa_add_graph_edge(ab, last_id, node_id, 1-last_new, add_read_id, read_id, read_ids_n);
                last_id = node_id; last_new = 0;
            }
        } else if (op == ABPOA_CINS || op == ABPOA_CSOFT_CLIP || op == ABPOA_CHARD_CLIP) {
            query_id = (abpoa_cigar[i] >> 34) & 0x3fffffff;
            len = (abpoa_cigar[i] >> 4) & 0x3fffffff;
            for (j = len-1; j >= 0; --j) { // XXX use dynamic id, instead of static query_id
                new_id = abpoa_add_graph_node(ab, seq[query_id-j]);
                abpoa_add_graph_edge(ab, last_id, new_id, 0, add_read_id, read_id, read_ids_n);
                last_id = new_id; last_new = 1;
            }
        } else if (op == ABPOA_CDEL) {
            // nothing;
            continue;
        }
    } 
    abpoa_add_graph_edge(ab, last_id, ABPOA_SINK_NODE_ID, 1-last_new, add_read_id, read_id, read_ids_n);
    graph->is_called_cons = graph->is_topological_sorted = 0;
    // abpoa_topological_sort(graph, abpt);
    return 0;
}

// reset allocated memery everytime init the graph
// * node
// * index_to_node_id/node_id_to_index/node_id_to_max_remain, max_pos_left/right
void abpoa_reset_graph(abpoa_t *ab, abpoa_para_t *abpt, int qlen) {
    int i, k, node_m;
    ab->abg->cons_l = 0; ab->abg->is_topological_sorted = ab->abg->is_called_cons = 0;
    for (i = 0; i < ab->abg->node_n; ++i) {
        ab->abg->node[i].in_edge_n = ab->abg->node[i].out_edge_n = ab->abg->node[i].aligned_node_n = 0;
        if (abpt->use_read_ids) {
            for (k = 0; k < ab->abg->node[i].read_ids_n; ++k) 
                ab->abg->node[i].read_ids[k] = 0;
        }
    }
    ab->abg->node_n = 2;
    if (qlen+2 > ab->abg->node_m) {
        node_m = qlen+2; kroundup32(node_m);
        ab->abg->node = (abpoa_node_t*)_err_realloc(ab->abg->node, node_m * sizeof(abpoa_node_t));
        for (i = ab->abg->node_m; i < node_m; ++i) 
            abpoa_set_graph_node(ab->abg, i);
        ab->abg->node_m = ab->abg->index_rank_m = node_m;
        ab->abg->index_to_node_id = (int*)_err_realloc(ab->abg->index_to_node_id, node_m * sizeof(int));
        ab->abg->node_id_to_index = (int*)_err_realloc(ab->abg->node_id_to_index, node_m * sizeof(int));
        if (abpt->out_msa || abpt->cons_agrm == ABPOA_HC || abpt->is_diploid) 
            ab->abg->node_id_to_msa_rank = (int*)_err_realloc(ab->abg->node_id_to_msa_rank, node_m * sizeof(int));
        if (abpt->wb >= 0) {
            ab->abg->node_id_to_max_pos_left = (int*)_err_realloc(ab->abg->node_id_to_max_pos_left, node_m * sizeof(int));
            ab->abg->node_id_to_max_pos_right = (int*)_err_realloc(ab->abg->node_id_to_max_pos_right, node_m * sizeof(int));
            ab->abg->node_id_to_max_remain = (int*)_err_realloc(ab->abg->node_id_to_max_remain, node_m * sizeof(int));
        } else if (abpt->zdrop > 0) {
            ab->abg->node_id_to_max_remain = (int*)_err_realloc(ab->abg->node_id_to_max_remain, node_m * sizeof(int));
        }
    }
    // fprintf(stderr, "qlen: %d, node_n: %d, node_m: %d\n", qlen, ab->abg->node_n, ab->abg->node_m);
}
