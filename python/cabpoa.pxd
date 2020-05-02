from libc.stdint cimport int8_t, uint8_t, int32_t, int64_t, uint32_t, uint64_t
from libc.stdio cimport FILE

cdef extern from "simd_instruction.h":
    int simd_check()

cdef extern from "abpoa.h":
    cdef int ABPOA_GLOBAL_MODE "ABPOA_GLOBAL_MODE"
    cdef int ABPOA_LOCAL_MODE "ABPOA_LOCAL_MODE"
    cdef int ABPOA_EXTEND_MODE "ABPOA_EXTEND_MODE"

    # gap mode
    cdef int ABPOA_LINEAR_GAP "ABPOA_LINEAR_GAP"
    cdef int ABPOA_AFFINE_GAP "ABPOA_AFFINE_GAP"
    cdef int ABPOA_CONVEX_GAP "ABPOA_CONVEX_GAP"

    cdef int ABPOA_EXTRA_B "ABPOA_EXTRA_B"
    cdef float ABPOA_EXTRA_F "ABPOA_EXTRA_F"

    cdef char *ABPOA_CIGAR_STR "ABPOA_CIGAR_STR"
    cdef int ABPOA_CMATCH "ABPOA_CMATCH"
    cdef int ABPOA_CINS "ABPOA_CINS"
    cdef int ABPOA_CDEL "ABPOA_CDEL"      
    cdef int ABPOA_CDIFF "ABPOA_CDIFF"     
    cdef int ABPOA_CSOFT_CLIP "ABPOA_CSOFT_CLIP"
    cdef int ABPOA_CHARD_CLIP "ABPOA_CHARD_CLIP"

    cdef int ABPOA_SRC_NODE_ID "ABPOA_SRC_NODE_ID"
    cdef int ABPOA_SINK_NODE_ID "ABPOA_SINK_NODE_ID"

    cdef int ABPOA_OUT_CONS "ABPOA_OUT_CONS"
    cdef int ABPOA_OUT_MSA "ABPOA_OUT_MSA"
    cdef int ABPOA_OUT_BOTH "ABPOA_OUT_BOTH"

    cdef int ABPOA_HB "ABPOA_HB"
    cdef int ABPOA_HC "ABPOA_HC"
    cdef int ABPOA_MF "ABPOA_MF"

    ctypedef struct abpoa_res_t:
        int n_cigar
        uint64_t *graph_cigar
        int node_s, node_e, query_s, query_e # for local and  extension mode
        int n_aln_bases, n_matched_bases


    ctypedef struct abpoa_para_t:
        int m
        int *mat # score matrix
        int match, mismatch, gap_open1, gap_open2, gap_ext1, gap_ext2
        int inf_min
        int wb # 1st part of extra band width
        float wf # 2nd part of extra band width. w=wb+wf*L (L is sequence length)
        int zdrop, end_bonus # from minimap2
        int simd_flag # available SIMD instruction
        # alignment mode
        uint8_t ret_cigar, rev_cigar, out_msa, out_cons, is_diploid, use_read_ids # mode: 0: global, 1: local, 2: extend
        char *out_pog
        int align_mode, gap_mode, cons_agrm
        double min_freq # for diploid data


    ctypedef struct abpoa_node_t:
        int node_id, index, rank
        int in_edge_n, in_edge_m
        int *in_id
        int out_edge_n, out_edge_m
        int *out_id
        int *out_weight
        uint64_t *read_ids
        int read_ids_n # for diploid
        int aligned_node_n, aligned_node_m
        int *aligned_node_id # mismatch; aligned node will have same rank
        int heaviest_weight, heaviest_out_id # for consensus
        uint8_t base # 0~m

    ctypedef struct abpoa_graph_t:
        abpoa_node_t *node
        int node_n, node_m, index_rank_m
        int *index_to_node_id
        int *node_id_to_index 
        int *node_id_to_min_rank
        int *node_id_to_max_rank
        int *node_id_to_max_remain
        int *node_id_to_msa_rank
        int cons_l, cons_m
        uint8_t *cons_seq
        uint8_t is_topological_sorted, is_called_cons, is_set_msa_rank

    ctypedef struct abpoa_simd_matrix_t:
        pass
    
    ctypedef struct abpoa_t:
        abpoa_graph_t *abg
        abpoa_simd_matrix_t *abm

    # init for abpoa parameters
    abpoa_para_t *abpoa_init_para()
    void abpoa_post_set_para(abpoa_para_t *abpt)
    void abpoa_free_para(abpoa_para_t *abpt)


    # init for alignment
    abpoa_t *abpoa_init()
    void abpoa_free(abpoa_t *ab, abpoa_para_t *abpt)

    # do msa for a set of input sequences
    int abpoa_msa(abpoa_t *ab, abpoa_para_t *abpt, int n_seqs, int *seq_lens, uint8_t **seqs, FILE *out_fp, uint8_t ***cons_seq, int **cons_l, int *cons_n, uint8_t ***msa_seq, int *msa_l)

    # clean alignment graph
    void abpoa_reset_graph(abpoa_t *ab, abpoa_para_t *abpt, int qlen)

    # align a sequence to a graph
    int abpoa_align_sequence_to_graph(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, abpoa_res_t *res)

    # add an alignment to a graph
    int abpoa_add_graph_node(abpoa_t *ab, uint8_t base)
    void abpoa_add_graph_edge(abpoa_t *ab, int from_id, int to_id, int check_edge, uint8_t add_read_id, int read_id, int read_ids_n)
    int abpoa_add_graph_alignment(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, int n_cigar, uint64_t *abpoa_cigar, int read_id, int tot_read_n)
    void abpoa_topological_sort(abpoa_t *ab, abpoa_para_t *abpt)

    # generate consensus sequence from graph
    # para:
    #   out_fp: consensus sequence output in FASTA format, set as NULL to disable
    #   cons_seq, cons_l, cons_n: store consensus sequences in variables, set cons_n as NULL to disable. 
    #     cons_seq: store consensus sequences
    #     cons_l: store consensus sequences length
    #     cons_n: store number of consensus sequences
    #     Note: cons_seq and cons_l need to be freed by user.
    int abpoa_generate_consensus(abpoa_t *ab, abpoa_para_t *abpt, int seq_n, FILE *out_fp, uint8_t ***cons_seq, int **cons_l, int *cons_n)
    # generate column multiple sequence alignment from graph
    void abpoa_generate_rc_msa(abpoa_t *ab, int seq_n, FILE *out_fp, uint8_t ***msa_seq, int *msa_l)

    # generate DOT graph plot 
    int abpoa_dump_pog(abpoa_t *ab, abpoa_para_t *abpt)
