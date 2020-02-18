#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "abpoa_align.h"
#include "utils.h"

#define TODO_M   0x1
#define TODO_X   0x2
#define TODO_MX  0x3
#define TODO_E   0x4
#define TODO_ME  0x5
#define TODO_F   0x8
#define TODO_MF  0x9
#define TODO_MEF 0xd

#define HAS_M    0x10
#define HAS_X    0x20
#define HAS_MX   0x30
#define HAS_E    0x40
#define HAS_ME   0x50
#define HAS_F    0x80


// xl, el, fl: current extension length of X, E and F
// h: max, e: vertical, f: horizontal
typedef struct {
    int64_t h:32, xl:32;
    int64_t e:32, el:32;
    int64_t f:32, fl:32; 
    uint8_t todo_map; // XXX need to stored seperately
} dp_matrix_t;

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
    abpt->ystop_iter = -1, abpt->ystop_min_processed_n = YSTOP_MIN_NUM, abpt->ystop_min_frac = YSTOP_MIN_FRAC;

    abpt->ret_cigar = 1;  // return cigar
    abpt->rev_cigar = 0;  // reverse cigar
    abpt->out_msa = 0;    // output msa
    abpt->out_cons = 0;   // output consensus sequence in msa
    abpt->multip = ABPOA_MULTIP; // muliple consensus for multiploid data
    abpt->min_fre = 0.3;
    abpt->cons_agrm = ABPOA_HB;   // consensus calling algorithm 
    abpt->use_read_ids = 0;
    abpt->out_pog= 0; // generate partial order graph

    // number of residue types
    abpt->m = 5; // nucleotides
    abpt->mat = NULL; // TODO score matrix for aa

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

void abpoa_free_para(abpoa_para_t *abpt) {
    if (abpt->mat != NULL) free(abpt->mat);
    free(abpt);
}
