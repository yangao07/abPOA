#ifndef POA_GRAPH_H
#define POA_GRAPH_H

#include <stdint.h>
#include "bpoa.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

//#define CIGAR_STR "MIDNSHP=XB"
//#define POA_GRAPH_CIGAR_STR "=XIDNSH"
//#define POA_GRAPH_CEQUAL 0
//#define POA_GRAPH_CMISMATCH 1
//#define POA_GRAPH_CINS 2
//#define POA_GRAPH_CDEL 3
//#define POA_GRAPH_CREF_SKIP 4
//#define POA_GRAPH_CCLIP 5

#define POA_SRC_NODE_ID 0
#define POA_SINK_NODE_ID 1

typedef struct {
    int node_id, index, rank;
    int in_edge_n, in_edge_m, *in_id;
    int out_edge_n, out_edge_m, *out_id, *out_weight;
    int aligned_node_n, aligned_node_m, *aligned_node_id; // mismatch; aligned node will have same rank
    int heavest_weight, heavest_out_id; // for consensus
    uint8_t base; // 0~m
    // ID, pos ???
} bpoa_node_t;

// TODO imitate bwt index, use uint64 for all variables
// XXX merge index_to_node_id and index_to_rank together uint64_t ???
typedef struct {
    bpoa_node_t *node; int node_n, node_m; 
    int *index_to_node_id, *index_to_min_rank, *index_to_max_rank, *index_to_min_remain, *index_to_max_remain;
    int *node_id_to_index, *node_id_to_min_rank, *node_id_to_max_rank, *node_id_to_min_remain, *node_id_to_max_remain;
    int min_rank_n, max_rank_n;
    int cons_l, cons_m; uint8_t *cons_seq;
    uint8_t is_topological_sorted:1, is_called_cons:1; 
} bpoa_graph_t;

// XXX max of in_edge is pow(2,30)
// for MATCH/MISMATCH: node_id << 34  | query_id << 4 | op
// for INSERTION:      query_id << 34 | op_len << 4   | op
// for DELETION:       node_id << 34  | op_len << 4   | op // op_len is always equal to 1
// for CLIP            query_id << 34 | op_len << 4   | op 
#define bpoa_cigar_t uint64_t 
                       
//typedef struct {
//    uint8_t op; // 0:match, 1:insertion, 2:deletion, 3:mismatch
//    int32_t len, m, *node_id, *query_id;
//} bpoa_graph_cigar_t; // alignment result of mapping a sequence to a graph


bpoa_node_t *bpoa_init_node(int n);
void bpoa_free_node(bpoa_node_t *node, int n);
bpoa_graph_t *bpoa_init_graph(int n);
void bpoa_free_graph(bpoa_graph_t *graph, int n);
int bpoa_align_sequence_with_graph(bpoa_graph_t *graph, uint8_t *query, int qlen, bpoa_para_t *bpt, int *n_cigar, bpoa_cigar_t **graph_cigar);
int bpoa_add_graph_alignment(bpoa_graph_t *graph, uint8_t *query, int qlen, int n_cigar, bpoa_cigar_t *bpoa_cigar, int *seq_node_ids, int *seq_node_ids_l);
int bpoa_generate_consensus(bpoa_graph_t *graph);
int bpoa_generate_multiple_sequence_alingment(bpoa_graph_t *graph, int **seq_node_ids, int *seq_node_ids_l, int seq_n, int output_consensu, FILE *fp);

static inline int bpoa_graph_node_id_to_index(bpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_index[node_id];
}

static inline int bpoa_graph_node_id_to_max_rank(bpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_max_rank[node_id];
}

static inline int bpoa_graph_node_id_to_min_rank(bpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_min_rank[node_id];
}

static inline int bpoa_graph_node_id_to_max_remain(bpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_max_remain[node_id];
}

static inline int bpoa_graph_node_id_to_min_remain(bpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_min_remain[node_id];
}

static inline int bpoa_graph_index_to_node_id(bpoa_graph_t *graph, int index_i) {
    if (index_i < 0 || index_i >= graph->node_n) err_fatal(__func__, "Wrong index: %d\n", index_i);
    return graph->index_to_node_id[index_i];
}

static inline int bpoa_graph_index_to_max_rank(bpoa_graph_t * graph, int index_i) {
    if (index_i < 0 || index_i >= graph->node_n) err_fatal(__func__, "Wrong index: %d\n", index_i);
    return graph->index_to_max_rank[index_i];
}

static inline int bpoa_graph_index_to_min_rank(bpoa_graph_t * graph, int index_i) {
    if (index_i < 0 || index_i >= graph->node_n) err_fatal(__func__, "Wrong index: %d\n", index_i);
    return graph->index_to_min_rank[index_i];
}

static inline int bpoa_graph_index_to_max_remain(bpoa_graph_t * graph, int index_i) {
    if (index_i < 0 || index_i >= graph->node_n) err_fatal(__func__, "Wrong index: %d\n", index_i);
    return graph->index_to_max_remain[index_i];
}

static inline int bpoa_graph_index_to_min_remain(bpoa_graph_t * graph, int index_i) {
    if (index_i < 0 || index_i >= graph->node_n) err_fatal(__func__, "Wrong index: %d\n", index_i);
    return graph->index_to_min_remain[index_i];
}

#ifdef __cplusplus
}
#endif

#endif
