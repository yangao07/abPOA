#include <stdio.h>
#include <stdlib.h>
#include "bpoa.h"
#include "bpoa_graph.h"
#include "bpoa_align.h"
#include "simd_bpoa_align.h"
#include "utils.h"
#include "kdq.h"


KDQ_INIT(int)
#define kdq_int_t kdq_t(int)

bpoa_node_t *bpoa_init_node(int n) {
    bpoa_node_t *node = (bpoa_node_t*)_err_calloc(n, sizeof(bpoa_node_t));
    return node;
}

void bpoa_set_graph_node(bpoa_graph_t *graph, int node_i) {
    graph->node[node_i].node_id = -1;
    graph->node[node_i].in_edge_n = 0;
    graph->node[node_i].in_edge_m = 0;
    graph->node[node_i].out_edge_n = 0;
    graph->node[node_i].out_edge_m = 0;
    graph->node[node_i].aligned_node_n = 0;
    graph->node[node_i].aligned_node_m = 0;

}

void bpoa_free_node(bpoa_node_t *node, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (node[i].in_edge_m > 0) free(node[i].in_id);
        if (node[i].out_edge_m > 0) { free(node[i].out_id); free(node[i].out_weight); }
        if (node[i].aligned_node_m > 0) free(node[i].aligned_node_id);
    }
    free(node);
}

// 0: in_edge, 1: out_edge
bpoa_graph_t *bpoa_realloc_graph_edge(bpoa_graph_t *graph, int io, int id) {
    if (io == 0) {
        _uni_realloc(graph->node[id].in_id, graph->node[id].in_edge_n, graph->node[id].in_edge_m, int);
    } else {
        int edge_m = graph->node[id].out_edge_m;
        _uni_realloc(graph->node[id].out_id, graph->node[id].out_edge_n, edge_m, int);
        _uni_realloc(graph->node[id].out_weight, graph->node[id].out_edge_n, graph->node[id].out_edge_m, int);
    }
    return graph;
}

bpoa_graph_t *bpoa_realloc_graph_node(bpoa_graph_t *graph) {
    if (graph->node_m <= 0) {
        graph->node_m = 1;
        graph->node = (bpoa_node_t*)_err_calloc(1, sizeof(bpoa_node_t));
    }
    if (graph->node_n == graph->node_m) {
        int i;
        graph->node_m <<= 1;
        graph->node = (bpoa_node_t*)_err_realloc(graph->node, graph->node_m * sizeof(bpoa_node_t));
        for (i = graph->node_m >> 1; i < graph->node_m; ++i) {
            bpoa_set_graph_node(graph, i);
        }
    }
    return graph;
}

bpoa_graph_t *bpoa_init_graph(int n) {
    bpoa_graph_t *graph = (bpoa_graph_t*)_err_calloc(n, sizeof(bpoa_graph_t));
    int i;
    for (i = 0; i < n; ++i) {
        graph[i].node_n = 2, graph[i].node_m = 2;
        graph[i].node = bpoa_init_node(2);
        graph[i].node[0].node_id = 0;
        graph[i].node[1].node_id = 1;
        graph[i].cons_l = graph[i].cons_m = 0;
        graph[i].is_topological_sorted = graph[i].is_called_cons = 0;
    }
    return graph;
}

void bpoa_free_graph(bpoa_graph_t *graph, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (graph[i].node_m > 0) bpoa_free_node(graph[i].node, graph[i].node_m);
        if (graph[i].cons_m > 0) free(graph[i].cons_seq);
        if (graph[i].node_n > 0) {
            free(graph[i].index_to_node_id);
            free(graph[i].index_to_min_rank);
            free(graph[i].index_to_max_rank);
            free(graph[i].index_to_min_remain);
            free(graph[i].index_to_max_remain);

            free(graph[i].node_id_to_index);
            free(graph[i].node_id_to_min_rank);
            free(graph[i].node_id_to_max_rank);
            free(graph[i].node_id_to_min_remain);
            free(graph[i].node_id_to_max_remain);
        }
    }
    free(graph);
}

// min_rank_to_index
// min_rank_to_node_id
// min_rank_index
// min_rank_n
/*
void bpoa_set_min_rank(bpoa_graph_t *graph) {
    int i, last_i = 0, min_rank;
    int *_min_rank_index = (int*)_err_calloc(graph->min_rank_n, sizeof(int));
    for (i = 0; i < graph->min_rank_n; ++i) {
        graph->min_rank_index[i] += last_i;
        _min_rank_index[i] = graph->min_rank_index[i];
        last_i = graph->min_rank_index[i];
    }
    for (i = 0; i < graph->node_n; ++i) {
        min_rank = graph->index_to_min_rank[i];
        graph->min_rank_to_index[_min_rank_index[min_rank]++] = i;
        min_rank = graph->node_id_to_min_rank[i];
        graph->min_rank_to_node_id[graph->min_rank_index[min_rank]++] = i;
    }
    free(_min_rank_index);
}*/

/*int bpoa_set_node_index_rank(bpoa_graph_t *graph, int src_id, int sink_id, int *in_degree) {
    int *id, cur_id, i, j, out_id, aligned_id, q_size, new_q_size;
    int rank = 0, rank_i = 0, index = 0;
    int out_min_rank = 0;
    int max_min_rank = 0;
    kdq_int_t *q = kdq_init_int();

    // node[q.id].in_degree equals 0
    // Breadth-First-Search
    kdq_push_int(q, src_id); q_size = 1; new_q_size = 0;
    graph->node_id_to_min_rank[src_id] = 0; 
    while (q_size > 0) {
        if ((id = kdq_shift_int(q)) == 0) err_fatal_simple("Error in queue.\n");
        cur_id = *id;
        graph->index_to_node_id[index] = cur_id;
        graph->index_to_rank[index] = rank; // XXX for aligned nodes, set same rank XXX
        graph->index_to_min_rank[index] = graph->node_id_to_min_rank[cur_id];
        if (graph->index_to_min_rank[index] > max_min_rank) max_min_rank = graph->index_to_min_rank[index];
        ++graph->min_rank_index[graph->index_to_min_rank[index]+1];
        graph->node_id_to_index[cur_id] = index++;
        graph->rank_to_node_id[rank_i++] = cur_id;
        graph->node_id_to_rank[cur_id] = rank;
        if (cur_id == sink_id) {
            graph->rank_n = rank+1;
            graph->min_rank_n = max_min_rank+1;
            bpoa_set_min_rank(graph);
            kdq_destroy_int(q);
            return 0;
        }

        out_min_rank = graph->node_id_to_min_rank[cur_id]+1;
        for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
            out_id = graph->node[cur_id].out_id[i];
            if (--in_degree[out_id] == 0) {
#ifdef __DEBUG__
                //printf("0-in: %d\n", out_id);
#endif
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
            // set min_rank
            if (graph->node_id_to_min_rank[out_id] > out_min_rank) {
                graph->node_id_to_min_rank[out_id] = out_min_rank;
            }
        }
        if (--q_size == 0) {
            graph->rank_index[rank] = rank_i;
            rank++;
            q_size = new_q_size;
            new_q_size = 0;
        }
    }
    return -1;
}*/

// TODO
// index_to_node_id, index_to_min_rank, index_max_rank, index_min_remain, index_to_max_remain
// node_id_to_index, node_id_to_min_rank, node_id_max_rank, node_id_min_remain, node_id_to_max_remain
int bpoa_set_node_index_rank(bpoa_graph_t *graph, int src_id, int sink_id, int *in_degree) {
    int *id, cur_id, i, j, out_id, in_id, aligned_id, q_size, new_q_size;
    int index = 0, out_min_rank, out_max_rank, in_min_remain, in_max_remain;
    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, src_id); q_size = 1; new_q_size = 0; // node[q.id].in_degree equals 0
    graph->node_id_to_min_rank[src_id] = graph->node_id_to_max_rank[src_id] = 0;

    while (q_size > 0) {
        if ((id = kdq_shift_int(q)) == 0) err_fatal_simple("Error in queue.\n");
        cur_id = *id;
        graph->index_to_node_id[index] = cur_id;
        // XXX for aligned nodes, set same rank XXX
        graph->index_to_min_rank[index] = graph->node_id_to_min_rank[cur_id];
        graph->index_to_max_rank[index] = graph->node_id_to_max_rank[cur_id];
        graph->node_id_to_index[cur_id] = index++;

        if (cur_id == sink_id) {
            kdq_destroy_int(q);
            goto set_remain;
        }

        out_min_rank = graph->node_id_to_min_rank[cur_id]+1;
        out_max_rank = graph->node_id_to_max_rank[cur_id]+1;

        for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
            out_id = graph->node[cur_id].out_id[i];
            if (--in_degree[out_id] == 0) {
#ifdef __DEBUG__
                //printf("0-in: %d\n", out_id);
#endif
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
              // set min_rank
              if (graph->node_id_to_min_rank[out_id] > out_min_rank) {
                  graph->node_id_to_min_rank[out_id] = out_min_rank;
              }
              // set max_rank
              if (graph->node_id_to_max_rank[out_id] < out_max_rank) {
                  graph->node_id_to_max_rank[out_id] = out_max_rank;
              }
        }
        if (--q_size == 0) {
            q_size = new_q_size;
            new_q_size = 0;
        }
    }
    return -1;

set_remain:
    graph->min_rank_n = graph->node_id_to_min_rank[sink_id]-1;
    graph->max_rank_n = graph->node_id_to_max_rank[sink_id]-1;
    
    graph->node_id_to_min_remain[sink_id] = graph->node_id_to_max_remain[sink_id] = -1;
    for (i = graph->node_n-1; i >= 0; --i) {
        cur_id = bpoa_graph_index_to_node_id(graph, i);
        graph->index_to_min_remain[i] = graph->node_id_to_min_remain[cur_id];
        graph->index_to_max_remain[i] = graph->node_id_to_max_remain[cur_id];
        in_min_remain = graph->node_id_to_min_remain[cur_id]+1;
        in_max_remain = graph->node_id_to_max_remain[cur_id]+1;
        for (j = 0; j < graph->node[cur_id].in_edge_n; ++j) {
            in_id = graph->node[cur_id].in_id[j];
            if (in_min_remain < graph->node_id_to_min_remain[in_id]) {
                graph->node_id_to_min_remain[in_id] = in_min_remain;
            }
            if (in_max_remain > graph->node_id_to_max_remain[in_id]) {
                graph->node_id_to_max_remain[in_id] = in_max_remain;
            }
        }
    }

#ifdef __DEBUG__
    for (i = 0; i < graph->node_n; ++i) {
        printf("%d, %d ==> %d\t%d\t%d\t%d\n", i, graph->index_to_node_id[i], graph->index_to_min_rank[i], graph->index_to_max_rank[i], graph->index_to_min_remain[i], graph->index_to_max_remain[i]);
    }
#endif
    return 0;
}

// XXX update instead of calculate from scratch
// 1. index_to_node_id
// 2. node_id_to_index
// 3. rank_to_node_id
// 4. node_id_to_rank
int bpoa_topological_sort(bpoa_graph_t *graph) {
    if (graph->node_n <= 0) {
        err_func_format_printf(__func__, "Empty graph.\n");
        return 0;
    }

    int i, node_n = graph->node_n;
    graph->index_to_node_id = (int*)_err_realloc(graph->index_to_node_id, node_n * sizeof(int));
    graph->index_to_min_rank = (int*)_err_realloc(graph->index_to_min_rank, node_n * sizeof(int));
    graph->index_to_max_rank = (int*)_err_realloc(graph->index_to_max_rank, node_n * sizeof(int));
    graph->index_to_min_remain = (int*)_err_realloc(graph->index_to_min_remain, node_n * sizeof(int));
    graph->index_to_max_remain = (int*)_err_realloc(graph->index_to_max_remain, node_n * sizeof(int));

    graph->node_id_to_index = (int*)_err_realloc(graph->node_id_to_index, node_n * sizeof(int));
    graph->node_id_to_min_rank = (int*)_err_realloc(graph->node_id_to_min_rank, node_n * sizeof(int));
    graph->node_id_to_max_rank = (int*)_err_realloc(graph->node_id_to_max_rank, node_n * sizeof(int));
    graph->node_id_to_min_remain = (int*)_err_realloc(graph->node_id_to_min_remain, node_n * sizeof(int));
    graph->node_id_to_max_remain = (int*)_err_realloc(graph->node_id_to_max_remain, node_n * sizeof(int));

    for (i = 0; i < node_n; ++i) {
        graph->node_id_to_min_rank[i] = graph->node_n;
        graph->node_id_to_max_rank[i] = 0;
        graph->node_id_to_min_remain[i] = graph->node_n;
        graph->node_id_to_max_remain[i] = 0;
    }

    int *in_degree = (int*)_err_malloc(graph->node_n * sizeof(int));
    for (i = 0; i < graph->node_n; ++i) in_degree[i] = graph->node[i].in_edge_n;

    // start from POA_SRC_NODE_ID to POA_SINK_NODE_ID
    if (bpoa_set_node_index_rank(graph, POA_SRC_NODE_ID, POA_SINK_NODE_ID, in_degree) < 0)
        err_fatal_simple("Failed to topological sort graph.\n");
    graph->is_topological_sorted = 1;
    free(in_degree);
    return 0;
}

// FIXME tend to include more bases
void bpoa_traverse_min_flow(bpoa_graph_t *graph,  int src_id, int sink_id, int *out_degree) {
    int *id, i, cur_id, in_id, out_id, max_id, max_w, max_out_i, min_w;
    kdq_int_t *q = kdq_init_int();
    kdq_push_int(q, sink_id);
    // reverse Breadth-First-Search
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == sink_id) {
            graph->node[sink_id].heavest_out_id = -1;
            graph->node[sink_id].heavest_weight = INF_32_MAX;
        } else {
            max_w = INF_32_MIN, max_id = -1, max_out_i = -1, min_w = INF_32_MAX;
            for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
                out_id = graph->node[cur_id].out_id[i];
                min_w = MIN_OF_TWO(graph->node[out_id].heavest_weight, graph->node[cur_id].out_weight[i]);
                if (max_w < min_w) {
                    max_w = min_w;
                    max_id = out_id;
                    max_out_i = i;
                } else if (max_w == min_w) {
                    if (graph->node[cur_id].out_weight[max_out_i] < graph->node[cur_id].out_weight[i]) {
                        // max_w = min_w;
                        max_id = out_id;
                        max_out_i = i;
                    }
                }
            }
            graph->node[cur_id].heavest_out_id = max_id;
            graph->node[cur_id].heavest_weight = max_w;
        }

        if (cur_id == src_id) {
            kdq_destroy_int(q);
            return;
        }
        for (i = 0; i < graph->node[cur_id].in_edge_n; ++i) {
            in_id = graph->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) 
                kdq_push_int(q, in_id);
        }
    }
}

void bpoa_generate_consensus_core(bpoa_graph_t *graph, int src_id, int sink_id) {
    bpoa_node_t *node;
    if (graph->max_rank_n >= graph->cons_m) {
        graph->cons_m = graph->max_rank_n + 1;
        graph->cons_seq = (uint8_t*)_err_realloc(graph->cons_seq, graph->cons_m * sizeof(uint8_t));
    }
    node = graph->node + graph->node[src_id].heavest_out_id;
    while (node->node_id != sink_id) {
#ifdef __DEBUG__
        //printf("%d: %d\n", node->node_id, node->heavest_weight);
#endif
        _uni_realloc(graph->cons_seq, graph->cons_l, graph->cons_m, uint8_t);
        graph->cons_seq[graph->cons_l++] = node->base;
        node = graph->node + node->heavest_out_id;
    }
}

// should always topological sort first, then generate consensus
int bpoa_generate_consensus(bpoa_graph_t *graph) {
    int i, *out_degree = (int*)_err_malloc(graph->node_n * sizeof(int));
    for (i = 0; i < graph->node_n; ++i) {
        out_degree[i] = graph->node[i].out_edge_n;
    }

    bpoa_traverse_min_flow(graph, POA_SRC_NODE_ID, POA_SINK_NODE_ID, out_degree);
    free(out_degree);
    // backtrack to generate consensu sequence
    bpoa_generate_consensus_core(graph, POA_SRC_NODE_ID, POA_SINK_NODE_ID);
    graph->is_called_cons = 1;
    return graph->cons_l;
}

int bpoa_generate_multiple_sequence_alingment(bpoa_graph_t *graph, int **seq_node_ids, int *seq_node_ids_l, int seq_n, int output_consensu, FILE *fp) {
    int i, j, k, cur_id, aligned_id, rank;
    char **msa_seq = (char**)_err_malloc((seq_n+output_consensu) * sizeof(char*));
    for (i = 0; i < seq_n+output_consensu; ++i) {
        msa_seq[i] = (char*)_err_malloc((graph->max_rank_n+1) * sizeof(char));
        for (j = 0; j < graph->max_rank_n; ++j) 
            msa_seq[i][j] = '-';
        msa_seq[i][j] = '\0';
    }
    for (i = 0; i < seq_n; ++i) {
        for (j = 0; j < seq_node_ids_l[i]; ++j) {
            cur_id = seq_node_ids[i][j];
            rank = bpoa_graph_node_id_to_max_rank(graph, cur_id);
            for (k = 0; k < graph->node[cur_id].aligned_node_n; ++k) {
                aligned_id = graph->node[cur_id].aligned_node_id[k];
                rank = MAX_OF_TWO(rank, bpoa_graph_node_id_to_max_rank(graph, aligned_id));
            }
            msa_seq[i][rank-1] = "ACGTN"[graph->node[cur_id].base];
        }
        fprintf(fp, "%s\n", msa_seq[i]);
    }

    if (output_consensu) {
        bpoa_node_t *node = graph->node + graph->node[POA_SRC_NODE_ID].heavest_out_id;
        while (node->node_id != POA_SINK_NODE_ID) {
            rank = bpoa_graph_node_id_to_max_rank(graph, node->node_id);
            for (k = 0; k < graph->node[node->node_id].aligned_node_n; ++k) {
                aligned_id = graph->node[node->node_id].aligned_node_id[k];
                rank = MAX_OF_TWO(rank, bpoa_graph_node_id_to_max_rank(graph, aligned_id));
            }
            msa_seq[seq_n][rank-1] = "ACGTN"[node->base];
            node = graph->node + node->heavest_out_id;
        }
        fprintf(fp, "%s\n", msa_seq[seq_n]);
    }

    for (i = 0; i < seq_n+output_consensu; ++i) free(msa_seq[i]); free(msa_seq);
    return 0;
}

int bpoa_get_aligned_id(bpoa_graph_t *graph, int node_id, uint8_t base) {
    int i;
    bpoa_node_t *node = graph->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        if (node[node[node_id].aligned_node_id[i]].base == base)
            return graph->node[node_id].aligned_node_id[i];
    }
    return -1;
}

void bpoa_add_graph_aligned_node1(bpoa_node_t *node, int aligned_id) {
    _uni_realloc(node->aligned_node_id, node->aligned_node_n, node->aligned_node_m, int);
    node->aligned_node_id[node->aligned_node_n++] = aligned_id;
}

void bpoa_add_graph_aligned_node(bpoa_graph_t *graph, int node_id, int aligned_id) {
    int i; bpoa_node_t *node = graph->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        bpoa_add_graph_aligned_node1(node + node[node_id].aligned_node_id[i], aligned_id);
        bpoa_add_graph_aligned_node1(node + aligned_id, node[node_id].aligned_node_id[i]);
    }
    bpoa_add_graph_aligned_node1(graph->node + node_id, aligned_id);
    bpoa_add_graph_aligned_node1(graph->node + aligned_id, node_id);
}

int bpoa_align_sequence_with_graph(bpoa_graph_t *graph, uint8_t *query, int qlen, bpoa_para_t *bpt, int *n_cigar, bpoa_cigar_t **graph_cigar) {
    if (graph->node_n <= 2 || qlen <= 0) { // empty graph or seq
        //err_func_format_printf(__func__, "graph node: %d\tquery: %d\n", graph->node_n, qlen);
        return -1;
    }
    int score = -1;
    if (bpt->align_mode == POA_GLOBAL_FLAG) {
        //score = bpoa_banded_global_align_sequence_with_graph(graph, query, qlen, bpt, n_cigar, graph_cigar);
        score = simd_bpoa_banded_global_align_sequence_with_graph(graph, query, qlen, bpt, n_cigar, graph_cigar);
    } else if (bpt->align_mode == POA_LOCAL_FLAG) {
        //score = bpoa_local_align_sequence_with_graph(graph, query, qlen, bpt, n_cigar, graph_cigar);
    } else if (bpt->align_mode == POA_EXTEND_FLAG) {
        //score = bpoa_extend_align_sequence_with_graph(graph, query, qlen, bpt, n_cigar, graph_cigar);
    } else {
        err_fatal(__func__, "Unknown align mode: %d\n", bpt->align_mode);
    }
    return score;
}

int bpoa_add_graph_node(bpoa_graph_t *graph, uint8_t base) {
    int node_id = graph->node_n;
    graph = bpoa_realloc_graph_node(graph);
    // add node
    graph->node[node_id].node_id = node_id;
    graph->node[node_id].base = base;

    ++graph->node_n;
    return node_id;
}

int bpoa_add_graph_edge(bpoa_graph_t *graph, int from_id, int to_id, int check_edge, int *seq_node_ids, int *seq_node_ids_l) {
    if (from_id < 0 || from_id >= graph->node_n || to_id < 0 || to_id >= graph->node_n) err_fatal(__func__, "node_n: %d\tfrom_id: %d\tto_id: %d\n", graph->node_n, from_id, to_id);

#ifdef __DEBUG__
    //printf("%d -> %d\n", from_id, to_id);
#endif
    // add to_id to seq_node_ids
    if (seq_node_ids && seq_node_ids_l) seq_node_ids[(*seq_node_ids_l)++] = to_id;

    if (check_edge) {
        int i;
        for (i = 0; i < graph->node[from_id].out_edge_n; ++i) {
            if (graph->node[from_id].out_id[i] == to_id) { // edge exists
                graph->node[from_id].out_weight[i]++; // update weight on existing edge
                // update label id
                return 0;
            }
        }
    }
    // add edge
    graph = bpoa_realloc_graph_edge(graph, 0, to_id); 
    graph->node[to_id].in_id[graph->node[to_id].in_edge_n] = from_id;
    ++graph->node[to_id].in_edge_n;
    graph = bpoa_realloc_graph_edge(graph, 1, from_id); 
    graph->node[from_id].out_id[graph->node[from_id].out_edge_n] = to_id;
    graph->node[from_id].out_weight[graph->node[from_id].out_edge_n] = 1; // initial weight for new edge
    ++graph->node[from_id].out_edge_n;
    return 0;
}

int bpoa_add_graph_sequence(bpoa_graph_t *graph, uint8_t *seq, int seq_l, int start, int end, int *seq_node_ids, int *seq_node_ids_l) {
    if (seq_l <= 0 || start >= seq_l || end <= start) err_fatal(__func__, "seq_l: %d\tstart: %d\tend: %d\n", seq_l, start, end);
    if (start < 0) start = 0; if (end > seq_l) end = seq_l;
    int node_id = bpoa_add_graph_node(graph, seq[start]);
    bpoa_add_graph_edge(graph, POA_SRC_NODE_ID, node_id, 0, seq_node_ids, seq_node_ids_l);
    int i; 
    for (i = start+1; i < end; ++i) {
        node_id = bpoa_add_graph_node(graph, seq[i]);
        bpoa_add_graph_edge(graph, node_id-1, node_id, 0, seq_node_ids, seq_node_ids_l);
    }
    bpoa_add_graph_edge(graph, node_id, POA_SINK_NODE_ID, 0, NULL, NULL);
    graph->is_topological_sorted = graph->is_called_cons = 0;
    bpoa_topological_sort(graph);
    return 0;
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
int bpoa_add_graph_alignment(bpoa_graph_t *graph, uint8_t *seq, int seq_l, int n_cigar, bpoa_cigar_t *bpoa_cigar, int *seq_node_ids, int *seq_node_ids_l) {
    if (graph->node_n == 2) { // empty graph
        bpoa_add_graph_sequence(graph, seq, seq_l, 0, seq_l, seq_node_ids, seq_node_ids_l);
        return 0;
    } else {
        if (graph->node_n < 2) {
            err_fatal(__func__, "Graph node: %d.\n", graph->node_n);
        } else if (n_cigar == 0) {
            err_fatal(__func__, "Empty graph cigar.\n");
        }
    }
    // normal graph, normal graph_cigar
    int i, j; int op, len, node_id, query_id, last_new = 0, last_id = POA_SRC_NODE_ID, new_id, aligned_id;
    for (i = 0; i < n_cigar; ++i) {
        op = bpoa_cigar[i] & 0xf;
        if (op == POA_CMATCH) {
            node_id = (bpoa_cigar[i] >> 34) & 0x3fffffff;
            query_id = (bpoa_cigar[i] >> 4) & 0x3fffffff;
            if (graph->node[node_id].base != seq[query_id]) { // mismatch
                // check if query base is identical to node_id's aligned node
                if ((aligned_id = bpoa_get_aligned_id(graph, node_id, seq[query_id])) >= 0) {
                    bpoa_add_graph_edge(graph, last_id, aligned_id, 1-last_new, seq_node_ids, seq_node_ids_l);
                    last_id = aligned_id; last_new = 0;
                } else {
                    new_id = bpoa_add_graph_node(graph, seq[query_id]);
                    bpoa_add_graph_edge(graph, last_id, new_id, 0, seq_node_ids, seq_node_ids_l);
                    last_id = new_id; last_new = 1;
                    // add new_id to node_id's aligned node
                    bpoa_add_graph_aligned_node(graph, node_id, new_id);
                }
            } else { // match
                bpoa_add_graph_edge(graph, last_id, node_id, 1-last_new, seq_node_ids, seq_node_ids_l);
                last_id = node_id; last_new = 0;
            }
        } else if (op == POA_CINS || op == POA_CSOFT_CLIP || op == POA_CHARD_CLIP) {
            query_id = (bpoa_cigar[i] >> 34) & 0x3fffffff;
            len = (bpoa_cigar[i] >> 4) & 0x3fffffff;
            for (j = len-1; j >= 0; --j) { // XXX use dynamic id, instead of static query_id
                new_id = bpoa_add_graph_node(graph, seq[query_id-j]);
                bpoa_add_graph_edge(graph, last_id, new_id, 0, seq_node_ids, seq_node_ids_l);
                last_id = new_id; last_new = 1;
            }
        } else if (op == POA_CDEL) {
            // nothing;
            continue;
        }
    } bpoa_add_graph_edge(graph, last_id, POA_SINK_NODE_ID, 1-last_new, NULL, NULL);
    graph->is_topological_sorted = graph->is_called_cons = 0;
    bpoa_topological_sort(graph);
    return 0;
}
