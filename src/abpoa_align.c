#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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
    abpt->bw = -1;        // disable bandwidth
    abpt->ystop_iter_n = -1, abpt->ystop_min_processed_n = YSTOP_MIN_NUM, abpt->ystop_min_frac = YSTOP_MIN_FRAC;

    abpt->ret_cigar = 1;  // return cigar
    abpt->rev_cigar = 0;  // reverse cigar
    abpt->out_msa = 0;    // output msa
    abpt->out_cons = 0;   // output consensus sequence in msa
    abpt->multip = ABPOA_MULTIP; // muliple consensus for multiploid data
    abpt->min_freq = 0.3;
    abpt->cons_agrm = ABPOA_HB;   // consensus calling algorithm 
    abpt->use_read_ids = 0;
    abpt->out_pog= NULL; // generate partial order graph

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
    if (abpt->cons_agrm == 2 || abpt->out_msa || abpt->multip > 1) {
        abpt->use_read_ids = 1;
        set_65536_table();
        if (abpt->cons_agrm == 2 || abpt->multip > 1) set_bit_table16();
    }
}

void abpoa_free_para(abpoa_para_t *abpt) {
    if (abpt->mat != NULL) free(abpt->mat);
    free(abpt);
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
    int i, only_cons = !abpt->out_msa && abpt->out_cons, check_ystop, ystop, ystop_n = 0, tot_n = n_seqs;
    abpoa_res_t res;
    for (i = 0; i < n_seqs; ++i) {
        res.graph_cigar = 0, res.n_cigar = 0;
        abpoa_align_sequence_to_graph(ab, abpt, seqs[i], seq_lens[i], &res);
        if (only_cons && abpt->ystop_iter_n > 0 && i >= abpt->ystop_min_processed_n) {
            abpt->ystop_min_wei = (i+1) * abpt->ystop_min_frac;
            check_ystop = 1;
        } else check_ystop = 0;
        ystop = abpoa_add_graph_alignment(ab, abpt, seqs[i], seq_lens[i], res.n_cigar, res.graph_cigar, check_ystop, i, n_seqs);
        if (res.n_cigar) free(res.graph_cigar);
        if (check_ystop) {
            if (ystop) ++ystop_n;
            else ystop_n = 0;
            if (ystop_n >= abpt->ystop_iter_n) {
                tot_n = i+1;
                break;
            }
        }
    }
    if (abpt->out_cons)
        abpoa_generate_consensus(ab, abpt->cons_agrm, abpt->multip, abpt->min_freq, tot_n, out_fp, cons_seq, cons_l, cons_n);
    if (abpt->out_msa) {
        abpoa_generate_rc_msa(ab, tot_n, out_fp, msa_seq, msa_l);
    }
    return 1;
}


