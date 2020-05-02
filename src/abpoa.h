#ifndef ABPOA_H
#define ABPOA_H

#include <stdint.h>
#include "simd_instruction.h"

char LogTable65536[65536];
char bit_table16[65536];

#define ABPOA_GLOBAL_MODE 0
#define ABPOA_LOCAL_MODE  1
#define ABPOA_EXTEND_MODE 2
//#define ABPOA_SEMI_MODE 3

// gap mode
#define ABPOA_LINEAR_GAP 0
#define ABPOA_AFFINE_GAP 1
#define ABPOA_CONVEX_GAP 2

#define ABPOA_EXTRA_B 10
#define ABPOA_EXTRA_F 0.01

#define ABPOA_CIGAR_STR "MIDXSH"
#define ABPOA_CMATCH     0
#define ABPOA_CINS       1
#define ABPOA_CDEL       2
#define ABPOA_CDIFF      3
#define ABPOA_CSOFT_CLIP 4
#define ABPOA_CHARD_CLIP 5

#define ABPOA_SRC_NODE_ID 0
#define ABPOA_SINK_NODE_ID 1

#define ABPOA_OUT_CONS 0
#define ABPOA_OUT_MSA 1
#define ABPOA_OUT_BOTH 2

#define ABPOA_HB 0
#define ABPOA_HC 1
#define ABPOA_MF 2

// NOTE: upper boundary of in_edge_n is pow(2,30)
// for MATCH/MISMATCH: node_id << 34  | query_id << 4 | op
// for INSERTION:      query_id << 34 | op_len << 4   | op
// for DELETION:       node_id << 34  | op_len << 4   | op // op_len is always equal to 1
// for CLIP            query_id << 34 | op_len << 4   | op 
#define abpoa_cigar_t uint64_t 

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int n_cigar; abpoa_cigar_t *graph_cigar;
    int node_s, node_e, query_s, query_e; // for local and  extension mode
    int n_aln_bases, n_matched_bases;
} abpoa_res_t;

typedef struct {
    int m; int *mat; // score matrix
    int match, mismatch, gap_open1, gap_open2, gap_ext1, gap_ext2; int inf_min;
    int wb; float wf; // extra band width
    int zdrop, end_bonus; // from minimap2
    int simd_flag; // available SIMD instruction
    // alignment mode
    uint8_t ret_cigar:1, rev_cigar:1, out_msa:1, out_cons:1, is_diploid:1, use_read_ids:1;
    char *out_pog;
    int align_mode, gap_mode, cons_agrm;
    double min_freq; // for multiploid data
} abpoa_para_t;

typedef struct {
    int node_id;
    int in_edge_n, in_edge_m, *in_id;
    int out_edge_n, out_edge_m, *out_id, *out_weight;
    uint64_t *read_ids; int read_ids_n; // for multiploid

    int aligned_node_n, aligned_node_m, *aligned_node_id; // mismatch; aligned node will have same rank
    // int heaviest_weight, heaviest_out_id; // for consensus
    uint8_t base; // 0~m
    // ID, pos ???
} abpoa_node_t;

typedef struct {
    abpoa_node_t *node; int node_n, node_m, index_rank_m; 
    int *index_to_node_id;
    int *node_id_to_index, *node_id_to_max_pos_left, *node_id_to_max_pos_right, *node_id_to_max_remain, *node_id_to_msa_rank;
    int cons_l, cons_m; uint8_t *cons_seq;
    uint8_t is_topological_sorted:1, is_called_cons:1, is_set_msa_rank:1;
} abpoa_graph_t;

typedef struct {
    SIMDi *s_mem; uint32_t s_msize; // qp, DP_HE, dp_f OR qp, DP_H, dp_f : based on (qlen, num_of_value, m, node_n)
    int *dp_beg, *dp_end, *dp_beg_sn, *dp_end_sn; int rang_m; // if band : based on (node_m)
} abpoa_simd_matrix_t;

typedef struct {
    abpoa_graph_t *abg;
    abpoa_simd_matrix_t *abm;
} abpoa_t;

// init for abpoa parameters
abpoa_para_t *abpoa_init_para(void);
void abpoa_post_set_para(abpoa_para_t *abpt);
void abpoa_free_para(abpoa_para_t *abpt);

// init for alignment
abpoa_t *abpoa_init(void);
void abpoa_free(abpoa_t *ab, abpoa_para_t *abpt);

// perform msa
int abpoa_msa(abpoa_t *ab, abpoa_para_t *abpt, int n_seqs, int *seq_lens, uint8_t **seqs, FILE *out_fp, uint8_t ***cons_seq, int **cons_l, int *cons_n, uint8_t ***msa_seq, int *msa_l);

// clean alignment graph
void abpoa_reset_graph(abpoa_t *ab, abpoa_para_t *abpt, int qlen);

// for development:
// align a sequence to a graph
int abpoa_align_sequence_to_graph(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, abpoa_res_t *res);

// add a node to a graph
// para:
//   base: 0123 for ACGT
int abpoa_add_graph_node(abpoa_t *ab, uint8_t base);

// add an edge to a graph
// para:
//   from_id/to_id: ids of from and to nodes
//   check_edge: set as 1 if this edge maybe alread exist and only need to update weight, set as 0 if the edge is new
//   add_read_id: set as 1 if read_id is used (to use row-column algorithm/generate MSA result/diploid consensus)
//   read_id: is of sequence
//   read_ids_n: size of read_id array, each one is 64-bit (1+(tot_read_n-1)/64)
int abpoa_add_graph_edge(abpoa_t *ab, int from_id, int to_id, int check_edge, uint8_t add_read_id, int read_id, int read_ids_n);

// add an alignment to a graph
// para:
//   query: 0123 for ACGT
//   qlen: query length
//   n_cigar/abpoa_cigar: from alignment result (abpoa_res_t)
//   read_id: id of sequence
//   tot_read_n: total number of sequence
int abpoa_add_graph_alignment(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, int n_cigar, abpoa_cigar_t *abpoa_cigar, int read_id, int tot_read_n);

// topological sortting of graph
void abpoa_topological_sort(abpoa_t *ab, abpoa_para_t *abpt);

// generate consensus sequence from graph
// para:
//   out_fp: consensus sequence output in FASTA format, set as NULL to disable
//   cons_seq, cons_l, cons_n: store consensus sequences in variables, set cons_n as NULL to disable. 
//     cons_seq: store consensus sequences
//     cons_l: store consensus sequences length
//     cons_n: store number of consensus sequences
//     Note: cons_seq and cons_l need to be freed by user.
int abpoa_generate_consensus(abpoa_t *ab, abpoa_para_t *abpt, int seq_n, FILE *out_fp, uint8_t ***cons_seq, int **cons_l, int *cons_n);

// generate column multiple sequence alignment from graph
void abpoa_generate_rc_msa(abpoa_t *ab, int seq_n, FILE *out_fp, uint8_t ***msa_seq, int *msa_l);

// generate DOT graph plot and dump graph into PDF/PNG format file
int abpoa_dump_pog(abpoa_t *ab, abpoa_para_t *abpt);

#ifdef __cplusplus
}
#endif

#endif
