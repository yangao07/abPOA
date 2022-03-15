/* example.c libabpoa usage example
   To compile: 
gcc -g example.c -I ./include -L ./lib -labpoa -lz -lm -o example
or:
gcc -g example.c -I ./include ./lib/libabpoa.a -lz -lm -o example
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "include/abpoa.h"

// for nt
// AaCcGgTtNn ==> 0,1,2,3,4
unsigned char nt4_table[256] = {
       0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
const char nt256_table[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

int main(void) {
    int i, j, n_seqs = 10;
    // char seqs[10][100] = {
    //     "CGTCAATCTATCGAAGCATACGCGGGCAGAGCCGAAGACCTCGGCAATCCA",
    //     "CCACGTCAATCTATCGAAGCATACGCGGCAGCCGAACTCGACCTCGGCAATCAC",
    //     "CGTCAATCTATCGAAGCATACGCGGCAGAGCCCGGAAGACCTCGGCAATCAC",
    //     "CGTCAATGCTAGTCGAAGCAGCTGCGGCAGAGCCGAAGACCTCGGCAATCAC",
    //     "CGTCAATCTATCGAAGCATTCTACGCGGCAGAGCCGACCTCGGCAATCAC",
    //     "CGTCAATCTAGAAGCATACGCGGCAAGAGCCGAAGACCTCGGCCAATCAC",
    //     "CGTCAATCTATCGGTAAAGCATACGCTCTGTAGCCGAAGACCTCGGCAATCAC",
    //     "CGTCAATCTATCTTCAAGCATACGCGGCAGAGCCGAAGACCTCGGCAATC",
    //     "CGTCAATGGATCGAGTACGCGGCAGAGCCGAAGACCTCGGCAATCAC",
    //     "CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT"
    //     };

    char seqs[10][100] = {
        "CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        "CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT"
        };

    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // alignment parameters
    // abpt->align_mode = 0; // 0:global 1:local, 2:extension
    // abpt->mat_fn = strdup("HOXD70.mtx"); abpt->use_score_matrix = 1; // score matrix instead of constant match/mismatch score
    // abpt->match = 2;      // match score
    // abpt->mismatch = 4;   // mismatch penalty
    // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
    // abpt->gap_open1 = 4;  // gap open penalty #1
    // abpt->gap_ext1 = 2;   // gap extension penalty #1
    // abpt->gap_open2 = 24; // gap open penalty #2
    // abpt->gap_ext2 = 1;   // gap extension penalty #2
                             // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    // abpt->bw = 10;        // extra band used in adaptive banded DP
    // abpt->bf = 0.01; 
     
    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
    abpt->progressive_poa = 1;
    abpt->max_n_cons = 2; // to generate 2 consensus sequences

    abpoa_post_set_para(abpt);

    // collect sequence length, trasform ACGT to 0123
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = strlen(seqs[i]);
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j)
            bseqs[i][j] = nt4_table[(int)seqs[i][j]];
    }

    // 1. directly output to stdout
    fprintf(stdout, "=== output to stdout ===\n");
    // perform abpoa-msa
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, stdout);

    // 2. output MSA alignment and  consensus sequence stored in (abpoa_cons_t *)
    abpoa_cons_t *abc = ab->abc;
    fprintf(stdout, "=== stored in variables ===\n");
    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (i = 0; i < abc->n_seq; ++i) {
        for (j = 0; j < abc->msa_len; ++j) {
            fprintf(stdout, "%c", nt256_table[abc->msa_base[i][j]]);
        }
        fprintf(stdout, "\n");
    }

    for (i = 0; i < abc->n_cons; ++i) {
        fprintf(stdout, ">Consensus_sequence");
        if (abc->n_cons > 1) {
            fprintf(stdout, "_%d ", i+1);
            for (j = 0; j < abc->clu_n_seq[i]; ++j) { // output read ids for each cluster/group
                fprintf(stdout, "%d", abc->clu_read_ids[i][j]);
                if (j != abc->clu_n_seq[i]-1) fprintf(stdout, ",");
            }
        }
        fprintf(stdout, "\n");
        for (j = 0; j < abc->cons_len[i]; ++j)
            fprintf(stdout, "%c", nt256_table[abc->cons_base[i][j]]);
        fprintf(stdout, "\n");
    }

    /* generate DOT partial order graph plot */
    abpt->out_pog = strdup("example.png"); // dump parital order graph to file
    if (abpt->out_pog != NULL) abpoa_dump_pog(ab, abpt);

    // free seq-related variables
    for (i = 0; i < n_seqs; ++i) free(bseqs[i]); free(bseqs); free(seq_lens);

    // free abpoa-related variables
    abpoa_free(ab); abpoa_free_para(abpt); 
    return 0;
}
