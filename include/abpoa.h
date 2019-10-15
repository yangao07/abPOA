#ifndef ABPOA_H
#define ABPOA_H

#include <stdint.h>
#include "simd_instruction.h"

char LogTable65536[65536];
char bit_table16[65536];
// XXX max of in_edge is pow(2,30)
// for MATCH/MISMATCH: node_id << 34  | query_id << 4 | op
// for INSERTION:      query_id << 34 | op_len << 4   | op
// for DELETION:       node_id << 34  | op_len << 4   | op // op_len is always equal to 1
// for CLIP            query_id << 34 | op_len << 4   | op 
#define abpoa_cigar_t uint64_t 

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int m; int *mat; // score matrix
    int match, mismatch, gap_open1, gap_open2, gap_ext1, gap_ext2; int inf_min;
    int bw; // band width
    int zdrop, end_bonus; // from minimap2
    int simd_flag; // available SIMD instruction
    // alignment mode
    uint8_t use_ada:1, ret_cigar:1, out_msa:1, out_cons:1, out_pog:1, use_read_ids:1; // mode: 0: global, 1: local, 2: extend
    int align_mode, gap_mode, cons_agrm;
    int multip; double min_fre; // for multiploid data
} abpoa_para_t;

typedef struct {
    int node_id, index, rank;
    int in_edge_n, in_edge_m, *in_id;
    int out_edge_n, out_edge_m, *out_id, *out_weight;
    uint64_t *read_ids; int read_ids_n; // for multiploid

    int aligned_node_n, aligned_node_m, *aligned_node_id; // mismatch; aligned node will have same rank
    int heaviest_weight, heaviest_out_id; // for consensus
    uint8_t base; // 0~m
    // ID, pos ???
} abpoa_node_t;

// TODO imitate bwt index, use uint64 for all variables
// XXX merge index_to_node_id and index_to_rank together uint64_t ???
typedef struct {
    abpoa_node_t *node; int node_n, node_m, index_rank_m; 
    int *index_to_node_id;
    int *node_id_to_index, *node_id_to_min_rank, *node_id_to_max_rank, *node_id_to_min_remain, *node_id_to_max_remain, *node_id_to_msa_rank;
    int cons_l, cons_m; uint8_t *cons_seq;
    uint8_t is_topological_sorted:1, is_called_cons:1, is_set_msa_rank:1;
} abpoa_graph_t;

typedef struct {
    SIMDi *s_mem; int s_msize; // qp, DP_HE, dp_f OR qp, DP_H, dp_f : based on (qlen, num_of_value, m, node_n)
    int *dp_beg, *dp_end, *dp_beg_sn, *dp_end_sn; int rang_m; // if band : based on (node_m)
    // int *pre_n, **pre_index; // pre_n, pre_index based on (node_n) TODO use in/out_id directly
} abpoa_simd_matrix_t;

typedef struct {
    abpoa_graph_t *abg;
    abpoa_simd_matrix_t *abm;
} abpoa_t;

// init for abpoa parameters
abpoa_para_t *abpoa_init_para(void);
void abpoa_free_para(abpoa_para_t *abpt);
void gen_simple_mat(int m, int *mat, int match, int mismatch);
void abpoa_set_gap_mode(abpoa_para_t *abpt);


// init for alignment
abpoa_t *abpoa_init(void);
void abpoa_free(abpoa_t *ab, abpoa_para_t *abpt);

void set_65536_table(void);
void set_bit_table16(void);

// clean alignment graph
void abpoa_reset_graph(abpoa_t *ab, int qlen, abpoa_para_t *abpt);

// align a sequence to a graph
int abpoa_align_sequence_with_graph(abpoa_t *ab, uint8_t *query, int qlen, abpoa_para_t *abpt, int *n_cigar, abpoa_cigar_t **graph_cigar);

// add an alignment to a graph
int abpoa_add_graph_alignment(abpoa_graph_t *graph, abpoa_para_t *abpt, uint8_t *query, int qlen, int n_cigar, abpoa_cigar_t *abpoa_cigar, int read_id, int read_ids_n);

// generate consensus sequence from graph
// para:
//   out_fp: consensus sequence output in FASTA format, set as NULL to disable
//   cons_seq, cons_l, cons_n: store consensus sequences in variables, set cons_n as NULL to disable. 
//     cons_seq: store consensus sequences
//     cons_l: store consensus sequences length
//     cons_n: store number of consensus sequences
//     Note: cons_seq and cons_l need to be freed by user.
int abpoa_generate_consensus(abpoa_graph_t *graph, uint8_t cons_agrm, int multip, double min_fre, int seq_n, FILE *out_fp, uint8_t **cons_seq, int *cons_l, int *cons_n);
// generate column multiple sequence alignment from graph
// store msa into msa[] // TODO
int abpoa_generate_multiple_sequence_alingment(abpoa_graph_t *graph, int seq_n, FILE *out_fp);

// generate DOT graph plot 
int abpoa_graph_visual(abpoa_graph_t *graph, char *dot_fn);
// int abpoa_main(const char *in_fn, int in_list, abpoa_para_t *abpt);

#ifdef __cplusplus
}
#endif

#endif
