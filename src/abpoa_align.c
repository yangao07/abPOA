#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "simd_abpoa_align.h"
#include "abpoa_align.h"
#include "utils.h"

void gen_simple_mat(int m, int *mat, int match, int mismatch) {
    int i, j;
    match = match < 0 ? -match : match;
    mismatch = mismatch > 0? -mismatch : mismatch;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j ? match : mismatch;
        mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = 0;
}

void abpoa_set_gap_mode(abpoa_para_t *abpt) {
    if (abpt->gap_open1 == 0) abpt->gap_mode = ABPOA_LINEAR_GAP;
    else if (abpt->gap_open1 > 0 && abpt->gap_open2 == 0) abpt->gap_mode = ABPOA_AFFINE_GAP;
    else abpt->gap_mode = ABPOA_CONVEX_GAP;
}

abpoa_para_t *abpoa_init_para(void) {
    abpoa_para_t *abpt = (abpoa_para_t*)_err_malloc(sizeof(abpoa_para_t));
    abpt->align_mode = ABPOA_GLOBAL_MODE;
    abpt->gap_mode = ABPOA_CONVEX_GAP;
    abpt->zdrop = -1;     // disable zdrop
    abpt->end_bonus = -1; // disable end bouns
    abpt->wb = ABPOA_EXTRA_B; // extra bandwidth
    abpt->wf = ABPOA_EXTRA_F; // extra bandwidth

    abpt->ret_cigar = 1;  // return cigar
    abpt->rev_cigar = 0;  // reverse cigar
    abpt->out_cons = 1;   // output consensus sequence in msa
    abpt->out_msa = 0;    // output msa
    abpt->is_diploid = 0; // diploid data
    abpt->min_freq = DIPLOID_MIN_FREQ; 
    abpt->cons_agrm = ABPOA_HB;   // consensus calling algorithm 
    abpt->use_read_ids = 0;
    abpt->out_pog = NULL; // dump partial order graph to file

    // number of residue types
    abpt->m = 5; // nucleotides TODO score matrix for aa
    abpt->mat = (int*)_err_malloc(abpt->m * abpt->m * sizeof(int));

    // score matrix
    abpt->match = ABPOA_MATCH;
    abpt->mismatch = ABPOA_MISMATCH;
    abpt->gap_open1 = ABPOA_GAP_OPEN1;
    abpt->gap_open2 = ABPOA_GAP_OPEN2;
    abpt->gap_ext1 = ABPOA_GAP_EXT1;
    abpt->gap_ext2 = ABPOA_GAP_EXT2;

    abpt->simd_flag = simd_check();

    return abpt;
}

void abpoa_post_set_para(abpoa_para_t *abpt) {
    gen_simple_mat(abpt->m, abpt->mat, abpt->match, abpt->mismatch);
    abpoa_set_gap_mode(abpt);
    if (abpt->cons_agrm == ABPOA_HC || abpt->out_msa || abpt->is_diploid) {
        abpt->use_read_ids = 1;
        set_65536_table();
        if (abpt->cons_agrm == ABPOA_HC || abpt->is_diploid) set_bit_table16();
    }
    if (abpt->align_mode == ABPOA_LOCAL_MODE) abpt->wb = -1;
}

void abpoa_free_para(abpoa_para_t *abpt) {
    if (abpt->mat != NULL) free(abpt->mat);
    if (abpt->out_pog != NULL) free(abpt->out_pog);
    free(abpt);
}

int abpoa_align_sequence_to_graph(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, abpoa_res_t *res) {
    if (ab->abg->node_n <= 2 || qlen <= 0) return -1;
    if (ab->abg->is_topological_sorted == 0) abpoa_topological_sort(ab, abpt);
    simd_abpoa_align_sequence_to_graph(ab, abpt, query, qlen, res);
    return 0;
}

// do msa for a set of input sequences
// @function: 
//    generate consensus sequence
//    generate rc-msa (row column multiple sequence alignment)
// @para:
//    ab/abpt: abpoa related variable and parameter 
//    seq_n: number of input sequences
//    seq_len: array of input sequence length, size: seq_n
//    seqs: array of input sequences, 0123 for ACGT, size: seq_n * seq_len[]
int abpoa_msa(abpoa_t *ab, abpoa_para_t *abpt, int n_seqs, int *seq_lens, uint8_t **seqs, FILE *out_fp, uint8_t ***cons_seq, int **cons_l, int *cons_n, uint8_t ***msa_seq, int *msa_l) {
    if ((!abpt->out_msa && !abpt->out_cons) || n_seqs <= 0) return 0;
    int i, tot_n = n_seqs;
    abpoa_reset_graph(ab, abpt, seq_lens[0]);
    abpoa_res_t res;
    for (i = 0; i < n_seqs; ++i) {
        res.graph_cigar = 0, res.n_cigar = 0;
        abpoa_align_sequence_to_graph(ab, abpt, seqs[i], seq_lens[i], &res);
        abpoa_add_graph_alignment(ab, abpt, seqs[i], seq_lens[i], res.n_cigar, res.graph_cigar, i, n_seqs);
        if (res.n_cigar) free(res.graph_cigar);
    }
    if (abpt->out_cons) {
        abpoa_generate_consensus(ab, abpt, tot_n, out_fp, cons_seq, cons_l, cons_n);
        if (ab->abg->is_called_cons == 0)
            err_printf("Warning: no consensus sequence generated.\n");
    }
    if (abpt->out_msa) {
        abpoa_generate_rc_msa(ab, tot_n, out_fp, msa_seq, msa_l);
    }
    return 1;
}
