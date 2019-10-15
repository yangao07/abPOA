#include <stdio.h>
#include <stdlib.h>
#include "abpoa_align.h"
#include "simd_abpoa_align.h"
#include "agglo_hier_clu.h"
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

void abpoa_free_node(abpoa_node_t *node, int n, abpoa_para_t *abpt) {
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
    graph->node_id_to_min_rank = NULL; graph->node_id_to_max_rank = NULL; graph->node_id_to_min_remain = NULL; graph->node_id_to_max_remain = NULL;
    return graph;
}

void abpoa_free_graph(abpoa_graph_t *graph, abpoa_para_t *abpt) {
    if (graph->node_m > 0) abpoa_free_node(graph->node, graph->node_m, abpt);
    if (graph->cons_m > 0) free(graph->cons_seq);

    if (graph->node_n > 0) {
        free(graph->index_to_node_id);
        free(graph->node_id_to_index);
        if (graph->node_id_to_msa_rank) free(graph->node_id_to_msa_rank);

        if (graph->node_id_to_min_rank) free(graph->node_id_to_min_rank);
        if (graph->node_id_to_max_rank) free(graph->node_id_to_max_rank);
        if (graph->node_id_to_min_remain) free(graph->node_id_to_min_remain);
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

int abpoa_BFS_set_node_index_rank(abpoa_graph_t *graph, int src_id, int sink_id, int *in_degree, int set_remain) {
    int *id, cur_id, i, j, out_id, in_id, aligned_id, q_size, new_q_size;
    int index = 0, in_min_remain, in_max_remain;
    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, src_id); q_size = 1; new_q_size = 0; // node[q.id].in_degree equals 0
    while (q_size > 0) {
        if ((id = kdq_shift_int(q)) == 0) err_fatal_simple("Error in queue.");
        cur_id = *id;
        graph->index_to_node_id[index] = cur_id;
        graph->node_id_to_index[cur_id] = index++;

        if (cur_id == sink_id) {
            kdq_destroy_int(q);
            if (set_remain)
                goto set_remain;
            else return 0;
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
    return -1;
set_remain:
    graph->node_id_to_min_remain[sink_id] = graph->node_id_to_max_remain[sink_id] = -1;
    for (i = graph->node_n-1; i >= 0; --i) {
        cur_id = abpoa_graph_index_to_node_id(graph, i);
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
    return 0;
}

// 1. index_to_node_id
// 2. node_id_to_index
// 3. node_id_to_rank
int abpoa_topological_sort(abpoa_graph_t *graph, abpoa_para_t *abpt) {
    if (graph->node_n <= 0) {
        err_func_format_printf(__func__, "Empty graph.\n");
        return 0;
    }
    int i, node_n = graph->node_n;
    if (node_n > graph->index_rank_m) {
        graph->index_rank_m = MAX_OF_TWO(node_n, graph->index_rank_m << 1);
        graph->index_to_node_id = (int*)_err_realloc(graph->index_to_node_id, graph->index_rank_m * sizeof(int));
        graph->node_id_to_index = (int*)_err_realloc(graph->node_id_to_index, graph->index_rank_m * sizeof(int));
        if (abpt->out_msa || abpt->cons_agrm == ABPOA_RC) graph->node_id_to_msa_rank = (int*)_err_realloc(graph->node_id_to_msa_rank, graph->index_rank_m * sizeof(int));
        if (abpt->bw >= 0) {
            graph->node_id_to_min_rank = (int*)_err_realloc(graph->node_id_to_min_rank, graph->index_rank_m * sizeof(int));
            graph->node_id_to_max_rank = (int*)_err_realloc(graph->node_id_to_max_rank, graph->index_rank_m * sizeof(int));
            graph->node_id_to_min_remain = (int*)_err_realloc(graph->node_id_to_min_remain, graph->index_rank_m * sizeof(int));
            graph->node_id_to_max_remain = (int*)_err_realloc(graph->node_id_to_max_remain, graph->index_rank_m * sizeof(int));
        }
    }
    if (abpt->bw >= 0) {
        for (i = 0; i < node_n; ++i) {
            graph->node_id_to_min_remain[i] = graph->node_n;
            graph->node_id_to_max_remain[i] = 0;
        }
    }

    int *in_degree = (int*)_err_malloc(graph->node_n * sizeof(int));
    for (i = 0; i < graph->node_n; ++i) in_degree[i] = graph->node[i].in_edge_n;

    // start from ABPOA_SRC_NODE_ID to ABPOA_SINK_NODE_ID
    if (abpoa_BFS_set_node_index_rank(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, in_degree, abpt->bw>=0) < 0)
        err_fatal_simple("Failed to topological sort graph.");
    graph->is_topological_sorted = 1;
    free(in_degree);
    return 0;
}

/* void abpoa_alloc_multip_path(abpoa_graph_t *graph, int multip, int read_ids_n) {
    int i, j; abpoa_node_t *node;
    for (i = 0; i < graph->node_n; ++i) {
        node = graph->node + i;
        node->path_n = 0;
        node->path_read_ids = (uint64_t**)_err_malloc(multip * sizeof(uint64_t*));
        for (j = 0; j < multip; ++j)
            node->path_read_ids[j] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
        node->path_out_i = (int*)_err_malloc(multip * sizeof(int));
        node->path_out_path_i = (int*)_err_malloc(multip * sizeof(int));
    }
}*/

/*void abpoa_free_multip_path(abpoa_graph_t *graph, int multip) {
    int i, j; abpoa_node_t *node;
    for (i = 0; i < graph->node_n; ++i) {
        node = graph->node + i;
        for (j = 0; j < multip; ++j) free(node->path_read_ids[j]);
        free(node->path_read_ids); free(node->path_out_i); free(node->path_out_path_i);
    }
}*/

/*void set_path_ids_with_edge(abpoa_node_t *node, int out_edge_i, int out_path_i, int multip, int read_ids_n) {
    int i, j; uint64_t b;
    if (node->path_n < multip) {
        for (i = 0; i < node->read_ids_n[out_edge_i]; ++i) {
            node->path_read_ids[node->path_n][i] = node->edge_read_ids[out_edge_i][i];
        }
        node->path_out_i[node->path_n] = out_edge_i; node->path_out_path_i[node->path_n] = out_path_i;
        ++node->path_n;
    } else {
        int min_w = INT32_MAX, min_i = -1, w, edge_w = node->out_weight[out_edge_i];
        for (i = 0; i < node->path_n; ++i) {
            w = 0;
            for (j = 0; j < read_ids_n; ++j) {
                b = node->path_read_ids[i][j];
                w += get_bit_cnt4(bit_table16, b);
            }
            if (w < min_w) {
                min_w = w;
                min_i = i;
            }
        }
        if (min_w < edge_w) {
            for (i = 0; i < node->read_ids_n[out_edge_i]; ++i)
                node->path_read_ids[min_i][i] = node->edge_read_ids[out_edge_i][i];
            node->path_out_i[min_i] = out_edge_i; node->path_out_path_i[min_i] = out_path_i;
        }
    }
}*/

/*void set_path_ids_with_path(abpoa_node_t *cur_node, int out_edge_i, abpoa_node_t *out_node, int out_path_i, int multip, int read_ids_n) {
    int i, j, iw, uw; double T = 0.8; // XXX// intersect, union
    uint64_t *ib = (uint64_t*)_err_malloc(sizeof(uint64_t) * read_ids_n);
    uint64_t *ub = (uint64_t*)_err_malloc(sizeof(uint64_t) * read_ids_n);
    iw = uw = 0;
    for (i = 0; i < cur_node->read_ids_n[out_edge_i]; ++i) {
        ib[i] = cur_node->edge_read_ids[out_edge_i][i] & out_node->path_read_ids[out_path_i][i];
        ub[i] = cur_node->edge_read_ids[out_edge_i][i] | out_node->path_read_ids[out_path_i][i];
        iw += get_bit_cnt4(bit_table16, ib[i]); uw += get_bit_cnt4(bit_table16, ub[i]);
    }
    if (iw > 0) {
        if (cur_node->path_n < multip) {
            if (uw * T <= iw) { // do union
                for (i = 0; i < cur_node->read_ids_n[out_edge_i]; ++i)
                    cur_node->path_read_ids[cur_node->path_n][i] = ub[i];
            } else { // do intersection
                for (i = 0; i < cur_node->read_ids_n[out_edge_i]; ++i)
                    cur_node->path_read_ids[cur_node->path_n][i] = ib[i];
            }
            cur_node->path_out_i[cur_node->path_n] = out_edge_i;
            cur_node->path_out_path_i[cur_node->path_n] = out_path_i;
            ++cur_node->path_n;
        } else {
            int min_w = INT32_MAX, min_i, w, path_w; uint64_t b;
            uint64_t *path_b;
            if (uw * T <= iw) {
                path_w = uw; path_b = ub;
            } else {
                path_w = iw; path_b = ib;
            }
            for (i = 0; i < cur_node->path_n; ++i) {
                w = 0;
                for (j = 0; j < read_ids_n; ++j) {
                    b = cur_node->path_read_ids[i][j];
                    w += get_bit_cnt4(bit_table16, b);
                }
                if (w < min_w) {
                    min_w = w;
                    min_i = i;
                }
            }
            if (min_w < path_w) {
                for (i = 0; i < read_ids_n; ++i)
                    cur_node->path_read_ids[min_i][i] = path_b[i];
                cur_node->path_out_i[min_i] = out_edge_i;
                cur_node->path_out_path_i[min_i] = out_path_i;
            }
        }
    }
    free(ib); free(ub);
}*/

typedef struct {
    int v, p;
} dou_t;

int dou_cmp(const void *a, const void *b) { return (((dou_t*)b)->v - (((dou_t*)a)->v)); }

/*void abpoa_multip_consensus(abpoa_graph_t *graph, int src_id, int sink_id, int multip, int read_ids_n) {
    int i, j; uint64_t b; 
    dou_t w[multip];
    for (i = 0; i < graph->node[src_id].path_n; ++i) {
        w[i].v = 0; w[i].p = i;
        for (j = 0; j < read_ids_n; ++j) {
            b = graph->node[src_id].path_read_ids[i][j];
            w[i].v += get_bit_cnt4(bit_table16, b);
        }
    }
    qsort(w, graph->node[src_id].path_n, sizeof(dou_t), dou_cmp);
    abpoa_node_t *node; int path_i, out_edge_i, out_id, out_path_i;
    for (i = 0; i < multip && i < graph->node[src_id].path_n; ++i) {
        fprintf(stdout, ">Consensus_%d\n", i+1);
        path_i = w[i].p; 
        out_edge_i = graph->node[src_id].path_out_i[path_i];
        out_id = graph->node[src_id].out_id[out_edge_i];
        out_path_i = graph->node[src_id].path_out_path_i[path_i];
        node = graph->node + out_id;
        while (node->node_id != sink_id) {
            fprintf(stdout, "%c", "ACGTN"[node->base]);
            out_edge_i = node->path_out_i[out_path_i];
            out_id = node->out_id[out_edge_i];
            out_path_i = node->path_out_path_i[out_path_i];
            node = graph->node + out_id;
        }
        fprintf(stdout, "\n");
    }
}*/

/*void abpoa_multip_min_flow(abpoa_graph_t *graph, int src_id, int sink_id, int *out_degree, int seq_n, int multip) {
    int read_ids_n = (seq_n-1)/64+1;
    abpoa_alloc_multip_path(graph, multip, read_ids_n);
    int *id, i, j, cur_id, in_id, out_id; 
    // reverse Breadth-First-Search
    kdq_int_t *q = kdq_init_int(); kdq_push_int(q, sink_id);
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id != sink_id) {
            for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
                out_id = graph->node[cur_id].out_id[i];
                if (out_id == sink_id) {
                    // get ids from edge, pre_path_i: -1
                    set_path_ids_with_edge(graph->node + cur_id, i, -1, multip, read_ids_n);
                } else {
                    for (j = 0; j < graph->node[out_id].path_n; ++j) {
                        // get ids from path, pre_path_i
                        set_path_ids_with_path(graph->node + cur_id, i,  graph->node + out_id, j, multip, read_ids_n);
                    }
                }
            }
        }
        if (cur_id == src_id) { kdq_destroy_int(q); break; }
        for (i = 0; i < graph->node[cur_id].in_edge_n; ++i) {
            in_id = graph->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) kdq_push_int(q, in_id);
        }
    }
    abpoa_multip_consensus(graph, src_id, sink_id, multip, read_ids_n);
    abpoa_free_multip_path(graph, multip);
}*/

int abpoa_DFS_set_msa_rank(abpoa_graph_t *graph, int src_id, int sink_id, int *in_degree) {
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
            return 0;
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
    graph->is_set_msa_rank = 1;
    return -1;
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

void abpoa_set_row_column_ids_weight(abpoa_graph_t *graph, uint64_t ***read_ids, int **rc_weight, int seq_n, int read_ids_n) {
    int i, j, k, rank, aligned_id;
    uint64_t b, one = 1, *whole_read_ids = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    for (i = 0; i < seq_n; ++i) {
        j = i / 64; b = i & 0x3f;
        whole_read_ids[j] |= (one << b);
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
            read_ids[rank-1][4][k] = whole_read_ids[k] ^ b;
        }
        rc_weight[rank-1][4] -= rc_weight[rank-1][graph->node[i].base];
    }
    free(whole_read_ids);
}

// param:
// cur_multip_n: cluster number in current multip_read_ids
// multip_n: cluster number of read_ids
// multip: max. clusters allowed
int add_multip_read_ids(int cur_multip_n, uint64_t **multip_read_ids, uint64_t **read_ids, int *multip_i, int multip_n, int read_ids_n) { // at most multip clusters
    int i, j;
    if (cur_multip_n == 0) {
        for (i = 0; i < multip_n; ++i) {
            for (j = 0; j < read_ids_n; ++j) multip_read_ids[i][j] = read_ids[multip_i[j]][j];
        }
        cur_multip_n = multip_n;
    } else {
        for (i = 0; i < multip_n; ++i) {
        }
        for (i = 0; i < read_ids_n; ++i) {
        }
    }
    return 0;
}

// select up to MULTIP optimal weighted nodes of each column
void abpoa_multip_row_column_ids_weight(int msa_l, uint64_t ***read_ids, int **rc_weight, int read_ids_n, int seq_n, int multip) {
    int i, j;
    uint64_t **multip_read_ids = (uint64_t**)_err_malloc(sizeof(uint64_t*) * sizeof(multip));
    for (i = 0; i < multip; ++i) multip_read_ids[i] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    int min_w = MAX_OF_TWO(1, seq_n / (multip + 1));
    int multip_n, w, cur_multip_n = 0;
    int *multip_i = (int*)_err_malloc(sizeof(int) * multip);
    for (i = 0; i < msa_l; ++i) {
        multip_n = 0;
        for (j = 0; j < 5; ++j) {
            w = rc_weight[i][j];
            if (w > min_w) {
                if (multip_n == multip) {
                    multip_n = 0;
                    break;
                }
                multip_i[++multip_n] = j;
            }
        }
        if (multip_n > 1 && multip_n <= multip) {
            add_multip_read_ids(cur_multip_n, multip_read_ids, read_ids[i], multip_i, multip_n, read_ids_n);
        }
    }
    for (i = 0; i < multip; ++i) free(multip_read_ids[i]); 
    free(multip_read_ids); free(multip_i);
}

void abpoa_set_weight_matrix(abpoa_graph_t *graph, int **weight_matrix, int seq_n) {
    int i = seq_n; i = graph->cons_l; weight_matrix[0][0] = 0; // TODO
}

void abpoa_multip_consensus(uint64_t ***read_ids, int **cluster_ids, int *cluster_ids_n, int multip, int read_ids_n, int msa_l, FILE *out_fp, uint8_t **_cons_seq, int *_cons_l, int *_cons_n) {
    int i, j, k, m, w, cnt, max_w, max_base, gap_w;
    uint64_t *read_ids_mask = (uint64_t*)_err_malloc(read_ids_n * sizeof(uint64_t)), one=1, b;
    uint8_t *cons_seq = (uint8_t*)_err_malloc(sizeof(uint8_t) * msa_l); int cons_l;
    int read_id, n;

    if (_cons_n) {
        *_cons_n = multip; 
        _cons_l = (int*)_err_malloc(sizeof(int) * multip);
        _cons_seq = (uint8_t **)_err_malloc(sizeof(uint8_t*) * multip);
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
        if (out_fp) for (j = 0; j < cons_l; ++j) fprintf(out_fp, "%c", "ACGT"[cons_seq[j]]); fprintf(out_fp, "\n");
        if (_cons_n) {
            _cons_l[i] = cons_l;
            _cons_seq[i] = (uint8_t*)_err_malloc(sizeof(uint8_t) * cons_l);
            for (j = 0; j < cons_l; ++j) _cons_seq[i][j] = cons_seq[j];
        }
    }
    free(read_ids_mask); free(cons_seq);
}

void abpoa_set_row_column_ids(abpoa_graph_t *graph, uint64_t ***read_ids, int seq_n, int read_ids_n) {
    int i, j, k, rank, aligned_id;
    uint64_t b, one = 1, *whole_read_ids = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    for (i = 0; i < seq_n; ++i) {
        j = i / 64; b = i & 0x3f;
        whole_read_ids[j] |= (one << b);
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
            read_ids[rank-1][graph->node[i].base][k] = b;
            read_ids[rank-1][4][k] = whole_read_ids[k] ^ b;
        }
    }
}

void abpoa_multip_cluster_row_column(abpoa_graph_t *graph, int src_id, int sink_id, int seq_n, int multip, double min_fre, FILE *out_fp, uint8_t **cons_seq, int *cons_l, int *cons_n) {
    abpoa_set_msa_rank(graph, src_id, sink_id);
    int i, j, msa_l = graph->node_id_to_msa_rank[sink_id]-1;
    int read_ids_n = (seq_n-1)/64+1;
    int **weight_matrix = (int**)_err_malloc(sizeof(int*) * seq_n);
    for (i = 0; i < seq_n; ++i) weight_matrix[i] = (int*)_err_malloc(sizeof(int) * seq_n);
    int **cluster_ids = (int**)_err_malloc(sizeof(int*) * multip);
    int *cluster_ids_n = (int*)_err_malloc(sizeof(int) * multip);
    for (i = 0; i < multip; ++i) cluster_ids[i] = (int*)_err_malloc(sizeof(int) * seq_n);
    uint64_t ***read_ids = (uint64_t***)_err_malloc(sizeof(uint64_t**) * msa_l);
    for (i = 0; i < msa_l; ++i) {
        read_ids[i] = (uint64_t**)_err_malloc(sizeof(uint64_t*) * 5);
        for (j = 0; j < 5; ++j) read_ids[i][j] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    }
    abpoa_set_weight_matrix(graph, weight_matrix, seq_n);
    agglo_hier_clu(weight_matrix, seq_n, multip, min_fre, cluster_ids, cluster_ids_n);
    abpoa_set_row_column_ids(graph, read_ids, seq_n, read_ids_n);
    abpoa_multip_consensus(read_ids, cluster_ids, cluster_ids_n, multip, read_ids_n, msa_l, out_fp, cons_seq, cons_l, cons_n);
    for (i = 0; i < seq_n; ++i) {
        free(weight_matrix[i]); 
    } free(weight_matrix);
    for (i = 0; i < multip; ++i) free(cluster_ids[i]);
    free(cluster_ids); free(cluster_ids_n);
    for (i = 0; i < msa_l; ++i) {
        for (j = 0; j < 5; ++j) {
            free(read_ids[i][j]);
        } free(read_ids[i]);
    } free(read_ids);
}

void output_consensus(abpoa_graph_t *graph, FILE *out_fp) {
    int i;
    fprintf(out_fp, ">Consensus_sequence\n");
    for (i = 0; i < graph->cons_l; ++i) {
        fprintf(out_fp, "%c", "ACGTN"[graph->cons_seq[i]]);
    } fprintf(out_fp, "\n");
}

void abpoa_store_consensus(abpoa_graph_t *graph, uint8_t **cons_seq, int *cons_l) {
    cons_seq = (uint8_t**)_err_malloc(sizeof(uint8_t*));
    cons_l = (int*)_err_malloc(sizeof(int));
    cons_seq[0] = (uint8_t*)_err_malloc(sizeof(uint8_t) * graph->cons_l);
    cons_l[0] = graph->cons_l;
    int i;
    for (i = 0; i < graph->cons_l; ++i) cons_seq[0][i] = graph->cons_seq[i];
}

void abpoa_generate_consensus_core(abpoa_graph_t *graph, int src_id, int sink_id) {
    abpoa_node_t *node;
    if (graph->node_n > graph->cons_m) {
        graph->cons_m = MAX_OF_TWO(graph->cons_m << 1, graph->node_n);
        graph->cons_seq = (uint8_t*)_err_realloc(graph->cons_seq, graph->cons_m * sizeof(uint8_t));
    }
    node = graph->node + graph->node[src_id].heaviest_out_id;
    while (node->node_id != sink_id) {
        _uni_realloc(graph->cons_seq, graph->cons_l, graph->cons_m, uint8_t);
        graph->cons_seq[graph->cons_l++] = node->base;
        node = graph->node + node->heaviest_out_id;
    }
}

void abpoa_traverse_min_flow(abpoa_graph_t *graph,  int src_id, int sink_id, int *out_degree) {
    int *id, i, cur_id, in_id, out_id, max_id, max_w, max_out_i, min_w;
    kdq_int_t *q = kdq_init_int(); kdq_push_int(q, sink_id);
    // reverse Breadth-First-Search
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == sink_id) {
            graph->node[sink_id].heaviest_out_id = -1;
            graph->node[sink_id].heaviest_weight = INT32_MAX;
        } else {
            max_w = INT32_MIN, max_id = -1, max_out_i = -1, min_w = INT32_MAX;
            for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
                out_id = graph->node[cur_id].out_id[i];
                min_w = MIN_OF_TWO(graph->node[out_id].heaviest_weight, graph->node[cur_id].out_weight[i]);
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
            graph->node[cur_id].heaviest_out_id = max_id;
            graph->node[cur_id].heaviest_weight = max_w;
        }

        if (cur_id == src_id) { 
            kdq_destroy_int(q); 
            goto MF_CONS; 
        }
        for (i = 0; i < graph->node[cur_id].in_edge_n; ++i) {
            in_id = graph->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) kdq_push_int(q, in_id);
        }
    }
MF_CONS:
    abpoa_generate_consensus_core(graph, src_id, sink_id);
}

// heaviest_bundling
// 1. argmax{cur->weight}
// 2. argmax{out_node->weight}
void abpoa_heaviest_bundling(abpoa_graph_t *graph,  int src_id, int sink_id, int *out_degree) {
    int *id, i, cur_id, in_id, out_id, max_id, max_w, out_w;
    kdq_int_t *q = kdq_init_int();
    kdq_push_int(q, sink_id);
    // reverse Breadth-First-Search
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == sink_id) {
            graph->node[sink_id].heaviest_out_id = -1;
            graph->node[sink_id].heaviest_weight = INT32_MAX;
        } else {
            max_w = INT32_MIN, max_id = -1;
            for (i = 0; i < graph->node[cur_id].out_edge_n; ++i) {
                out_id = graph->node[cur_id].out_id[i];
                out_w = graph->node[cur_id].out_weight[i];
                if (max_w < out_w) {
                    max_w = out_w;
                    max_id = out_id;
                } else if (max_w == out_w) {
                    if (graph->node[max_id].heaviest_weight < graph->node[out_id].heaviest_weight) {
                        max_id = out_id;
                    }
                }
            }
            graph->node[cur_id].heaviest_out_id = max_id;
            graph->node[cur_id].heaviest_weight = max_w;
        }
        if (cur_id == src_id) {
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
    abpoa_generate_consensus_core(graph, src_id, sink_id);
}

void abpoa_row_column_consensus(abpoa_graph_t *graph, int **rc_weight, int msa_l, int seq_n) {
    if (graph->cons_m < msa_l+1) {
        graph->cons_m = MAX_OF_TWO(graph->cons_m << 1, msa_l+1);
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


void abpoa_row_column(abpoa_graph_t *graph, int src_id, int sink_id, int seq_n) {
    abpoa_set_msa_rank(graph, src_id, sink_id);

    int i, msa_l = graph->node_id_to_msa_rank[sink_id] - 1;
    int **rc_weight = (int**)_err_malloc(sizeof(int*) * msa_l);
    for (i = 0; i < msa_l; ++i) {
        rc_weight[i] = (int*)_err_calloc(5, sizeof(int)); // ACGT
        rc_weight[i][4] = seq_n;
    } 
    abpoa_set_row_column_weight(graph, rc_weight);

    abpoa_row_column_consensus(graph, rc_weight, msa_l, seq_n);
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

int set_clu_read_ids(int **clu_read_ids, int *clu_read_ids_n, int seq_n, uint64_t ***read_ids, uint8_t **het_cons_base, int *het_pos, int het_n) {
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
        if (diff[i][1] < diff[i][0]) {
            clu_read_ids[1][clu_read_ids_n[1]++] = i;
        } else {
            clu_read_ids[0][clu_read_ids_n[0]++] = i;
        }
    }
    int clu_n = 0;
    for (i = 0; i < 2; ++i) {
        if (clu_read_ids_n[i] > 0) ++clu_n;
    }
    for (i = 0; i < seq_n; ++i) free(diff[i]); free(diff);
    return clu_n;
}

int abpoa_diploid_ids(uint64_t ***read_ids, int **rc_weight, int msa_l, int seq_n, double min_fre, int read_ids_n, int **clu_read_ids, int *clu_read_ids_n) {
    int i, j;
    uint64_t **het_read_ids = (uint64_t**)_err_malloc(sizeof(uint64_t*) * 2 );
    for (i = 0; i < 2; ++i) het_read_ids[i] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    int *het_pos = (int*)_err_malloc(sizeof(int) * msa_l), het_n = 0;
    uint8_t **het_cons_base = (uint8_t**)_err_malloc(sizeof(uint8_t*) * msa_l);
    for (i = 0; i < msa_l; ++i) het_cons_base[i] = (uint8_t*)_err_malloc(sizeof(uint8_t) * 2);
    int min_w = MAX_OF_TWO(1, seq_n * min_fre);
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
        clu_n = set_clu_read_ids(clu_read_ids, clu_read_ids_n, seq_n, read_ids, het_cons_base, het_pos, het_n);
    }
    
    for (i = 0; i < msa_l; ++i) free(het_cons_base[i]);
    for (i = 0; i < 2; ++i) free(het_read_ids[i]); 
    free(het_cons_base); free(het_read_ids); free(het_pos);
    return clu_n;
}

void abpoa_diploid_row_column(abpoa_graph_t *graph, int src_id, int sink_id, int seq_n, double min_fre, FILE *out_fp, uint8_t **cons_seq, int *cons_l, int *cons_n) {
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
    abpoa_set_row_column_ids_weight(graph, read_ids, rc_weight, seq_n, read_ids_n);

    int **clu_read_ids = (int**)_err_malloc(sizeof(int*) * 2), *clu_read_ids_n;
    for (i = 0; i < 2; ++i) clu_read_ids[i] = (int*)_err_malloc(sizeof(int) * seq_n);
    clu_read_ids_n = (int*)_err_calloc(2, sizeof(int));
    int clu_n = abpoa_diploid_ids(read_ids, rc_weight, msa_l, seq_n, min_fre, read_ids_n, clu_read_ids, clu_read_ids_n);
    if (clu_n == 1) {
        abpoa_row_column_consensus(graph, rc_weight, msa_l, seq_n);
        if (out_fp) output_consensus(graph, out_fp);
        if (cons_n) {
            *cons_n = 1; abpoa_store_consensus(graph, cons_seq, cons_l);
        }
    } else abpoa_multip_consensus(read_ids, clu_read_ids, clu_read_ids_n, clu_n, read_ids_n, msa_l, out_fp, cons_seq, cons_l, cons_n);

    for (i = 0; i < msa_l; ++i) {
        for (j = 0; j < 5; ++j) {
            free(read_ids[i][j]);
        } free(read_ids[i]); free(rc_weight[i]);
    } free(read_ids); free(rc_weight);
    for (i = 0; i < 2; ++i) free(clu_read_ids[i]); free(clu_read_ids); free(clu_read_ids_n);
}

// should always topological sort first, then generate consensus
int abpoa_generate_consensus(abpoa_graph_t *graph, uint8_t cons_agrm, int multip, double min_fre, int seq_n, FILE *out_fp, uint8_t **cons_seq, int *cons_l, int *cons_n) {
    int i, *out_degree = (int*)_err_malloc(graph->node_n * sizeof(int));
    for (i = 0; i < graph->node_n; ++i) {
        out_degree[i] = graph->node[i].out_edge_n;
    }

    if (multip == 1) {
        if (cons_agrm == ABPOA_HB) abpoa_heaviest_bundling(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree);
        else if (cons_agrm == ABPOA_MF) abpoa_traverse_min_flow(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree); 
        else if (cons_agrm == ABPOA_RC) abpoa_row_column(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, seq_n);
        else err_fatal(__func__, "Unknown consensus calling algorithm: %d.", cons_agrm);
        if (out_fp) output_consensus(graph, out_fp);
        if (cons_n) {
            *cons_n = 1; abpoa_store_consensus(graph, cons_seq, cons_l);
        }
    } else if (multip == 2) {
        abpoa_diploid_row_column(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, seq_n, min_fre, out_fp, cons_seq, cons_l, cons_n);
    } else if (multip > 2) {
        abpoa_multip_cluster_row_column(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, seq_n, multip, min_fre, out_fp, cons_seq, cons_l, cons_n);
    }
    graph->is_called_cons = 1; free(out_degree);
    return graph->cons_l;
}

void abpoa_set_msa_seq(abpoa_node_t node, int rank, char **msa_seq) {
    int i, b, read_id; char base = "ACGTN"[node.base];
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

int abpoa_generate_multiple_sequence_alingment(abpoa_graph_t *graph, int seq_n, FILE *out_fp) {
    abpoa_set_msa_rank(graph, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);

    int i, j, k, aligned_id, rank;
    char **msa_seq = (char**)_err_malloc(seq_n * sizeof(char*));
    int msa_l = graph->node_id_to_msa_rank[ABPOA_SINK_NODE_ID]-1;

    for (i = 0; i < seq_n; ++i) {
        msa_seq[i] = (char*)_err_malloc((msa_l+1) * sizeof(char));
        for (j = 0; j < msa_l; ++j) 
            msa_seq[i][j] = '-';
        msa_seq[i][j] = '\0';
    }

    fprintf(out_fp, ">Multiple_sequence_alignment\n");
    for (i = 2; i < graph->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(graph, i);
        for (k = 0; k < graph->node[i].aligned_node_n; ++k) {
            aligned_id = graph->node[i].aligned_node_id[k];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(graph, aligned_id));
        }
        // assign seq
        abpoa_set_msa_seq(graph->node[i], rank, msa_seq);
    }
    for (i = 0; i < seq_n; ++i) fprintf(out_fp, "%s\n", msa_seq[i]);

    for (i = 0; i < seq_n; ++i) free(msa_seq[i]); free(msa_seq);
    return 0;
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

void abpoa_add_graph_aligned_node(abpoa_graph_t *graph, int node_id, int aligned_id) {
    int i; abpoa_node_t *node = graph->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        abpoa_add_graph_aligned_node1(node + node[node_id].aligned_node_id[i], aligned_id);
        abpoa_add_graph_aligned_node1(node + aligned_id, node[node_id].aligned_node_id[i]);
    }
    abpoa_add_graph_aligned_node1(graph->node + node_id, aligned_id);
    abpoa_add_graph_aligned_node1(graph->node + aligned_id, node_id);
}

int abpoa_align_sequence_with_graph(abpoa_t *ab, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar) {
    if (ab->abg->node_n <= 2 || qlen <= 0) return -1;
    else return simd_abpoa_align_sequence_with_graph(ab, query, qlen, abpt, n_cigar, graph_cigar);
}

//TODO
void abpoa_set_read_id(uint64_t *read_ids, int read_id) {
    int n = read_id / 64;
    // if ((n = read_id / 64) >= read_ids_n) err_fatal("abpoa_set_read_id", "Unexpected read id: %d\n", read_id);
    uint64_t one = 1; int b = read_id & 0x3f;
    read_ids[n] |= (one << b);
}

int abpoa_add_graph_node(abpoa_graph_t *graph, uint8_t base) {
    int node_id = graph->node_n;
    graph = abpoa_realloc_graph_node(graph);
    // add node
    graph->node[node_id].base = base;
    ++graph->node_n;
    return node_id;
}

int abpoa_add_graph_edge(abpoa_graph_t *graph, int from_id, int to_id, int check_edge, uint8_t add_read_id, int read_id, int read_ids_n) {
    if (from_id < 0 || from_id >= graph->node_n || to_id < 0 || to_id >= graph->node_n) err_fatal(__func__, "node_n: %d\tfrom_id: %d\tto_id: %d.", graph->node_n, from_id, to_id);
    int out_edge_n = graph->node[from_id].out_edge_n;
    if (check_edge) {
        int i;
        for (i = 0; i < out_edge_n; ++i) {
            if (graph->node[from_id].out_id[i] == to_id) { // edge exists
                graph->node[from_id].out_weight[i]++; // update weight on existing edge
                // update label id
                goto ADD_READ_ID;
                return 0;
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
    return 0;
}

int abpoa_add_graph_sequence(abpoa_graph_t *graph, abpoa_para_t *abpt, uint8_t *seq, int seq_l, int start, int end, uint8_t add_read_id, int read_id, int read_ids_n) {
    if (seq_l <= 0 || start >= seq_l || end <= start) err_fatal(__func__, "seq_l: %d\tstart: %d\tend: %d.", seq_l, start, end);
    if (start < 0) start = 0; if (end > seq_l) end = seq_l;
    int node_id = abpoa_add_graph_node(graph, seq[start]);
    abpoa_add_graph_edge(graph, ABPOA_SRC_NODE_ID, node_id, 0, add_read_id, read_id, read_ids_n);
    int i; 
    for (i = start+1; i < end; ++i) {
        node_id = abpoa_add_graph_node(graph, seq[i]);
        abpoa_add_graph_edge(graph, node_id-1, node_id, 0, add_read_id, read_id, read_ids_n);
    }
    abpoa_add_graph_edge(graph, node_id, ABPOA_SINK_NODE_ID, 0, add_read_id, read_id, read_ids_n);
    graph->is_topological_sorted = graph->is_called_cons = graph->is_set_msa_rank = 0;
    abpoa_topological_sort(graph, abpt);
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
int abpoa_add_graph_alignment(abpoa_graph_t *graph, abpoa_para_t *abpt, uint8_t *seq, int seq_l, int n_cigar, abpoa_cigar_t *abpoa_cigar, int read_id, int read_ids_n) {
    uint8_t add_read_id = abpt->use_read_ids;
    if (graph->node_n == 2) { // empty graph
        abpoa_add_graph_sequence(graph, abpt, seq, seq_l, 0, seq_l, add_read_id, read_id, read_ids_n);
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
                    abpoa_add_graph_edge(graph, last_id, aligned_id, 1-last_new, add_read_id, read_id, read_ids_n);
                    last_id = aligned_id; last_new = 0;
                } else {
                    new_id = abpoa_add_graph_node(graph, seq[query_id]);
                    abpoa_add_graph_edge(graph, last_id, new_id, 0, add_read_id, read_id, read_ids_n);
                    last_id = new_id; last_new = 1;
                    // add new_id to node_id's aligned node
                    abpoa_add_graph_aligned_node(graph, node_id, new_id);
                }
            } else { // match
                abpoa_add_graph_edge(graph, last_id, node_id, 1-last_new, add_read_id, read_id, read_ids_n);
                last_id = node_id; last_new = 0;
            }
        } else if (op == ABPOA_CINS || op == ABPOA_CSOFT_CLIP || op == ABPOA_CHARD_CLIP) {
            query_id = (abpoa_cigar[i] >> 34) & 0x3fffffff;
            len = (abpoa_cigar[i] >> 4) & 0x3fffffff;
            for (j = len-1; j >= 0; --j) { // XXX use dynamic id, instead of static query_id
                new_id = abpoa_add_graph_node(graph, seq[query_id-j]);
                abpoa_add_graph_edge(graph, last_id, new_id, 0, add_read_id, read_id, read_ids_n);
                last_id = new_id; last_new = 1;
            }
        } else if (op == ABPOA_CDEL) {
            // nothing;
            continue;
        }
    } abpoa_add_graph_edge(graph, last_id, ABPOA_SINK_NODE_ID, 1-last_new, add_read_id, read_id, read_ids_n);
    graph->is_topological_sorted = graph->is_called_cons = 0;
    abpoa_topological_sort(graph, abpt);
    return 0;
}

// reset allocated memery everytime init the graph
// * node
// * index_to_node_id/node_id_to_index/node_id_to_max/min_rank/remain
void abpoa_reset_graph(abpoa_t *ab, int qlen, abpoa_para_t *abpt) {
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
    if ((qlen + 2) * 2 > ab->abg->node_m) {
        node_m = MAX_OF_TWO(ab->abg->node_m << 1, (qlen+2) * 2);
        ab->abg->node = (abpoa_node_t*)_err_realloc(ab->abg->node, node_m * sizeof(abpoa_node_t));
        for (i = ab->abg->node_m; i < node_m; ++i) 
            abpoa_set_graph_node(ab->abg, i);
        ab->abg->node_m = ab->abg->index_rank_m = node_m;
        ab->abg->index_to_node_id = (int*)_err_realloc(ab->abg->index_to_node_id, node_m * sizeof(int));
        ab->abg->node_id_to_index = (int*)_err_realloc(ab->abg->node_id_to_index, node_m * sizeof(int));
        ab->abg->node_id_to_msa_rank = (int*)_err_realloc(ab->abg->node_id_to_msa_rank, node_m * sizeof(int));
        if (abpt->bw >= 0) {
            ab->abg->node_id_to_min_rank = (int*)_err_realloc(ab->abg->node_id_to_min_rank, node_m * sizeof(int));
            ab->abg->node_id_to_max_rank = (int*)_err_realloc(ab->abg->node_id_to_max_rank, node_m * sizeof(int));
            ab->abg->node_id_to_min_remain = (int*)_err_realloc(ab->abg->node_id_to_min_remain, node_m * sizeof(int));
            ab->abg->node_id_to_max_remain = (int*)_err_realloc(ab->abg->node_id_to_max_remain, node_m * sizeof(int));
        }
    }
}
