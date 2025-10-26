/* incre_example.c: incremental MSA example using abPOA library
   To compile: 
gcc -g incre_example.c -I ./include -L ./lib -labpoa -lz -lm -o incre_example
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
    int i, j, n_seqs = 10, first_part = 5, second_part=3, third_part=2; // first_part: number of sequences in the first incremental MSA
                                           // append the remaining sequences later
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
    abpt->max_n_cons = 2; // to generate 2 consensus sequences

    abpoa_post_set_para(abpt);

    // collect sequence length, trasform ACGT to 0123
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = strlen(seqs[i]);
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j) {
            bseqs[i][j] = nt4_table[(int)seqs[i][j]];
        }
    }
    // abpt->verbose = ABPOA_DEBUG_VERBOSE;
    // 1. align the first part of sequences
    abpoa_msa(ab, abpt, first_part, NULL, seq_lens, bseqs, NULL, stdout);
    // 2. clean MSA & consensus-related memory to avoid leak, always call before incremental MSA
    abpoa_clean_msa_cons(ab);

    // 3. append the second part of sequences one by one
    printf("\nExisting graph has %d sequences, %d nodes. Now append the remaining %d sequences one by one.\n\n", ab->abs->n_seq, ab->abg->node_n, second_part);
    abpt->out_msa = 0, abpt->out_cons = 0; // disable output for intermediate append
    for (i = 0; i < second_part; ++i) {
        abpoa_msa(ab, abpt, 1, NULL, seq_lens + first_part + i, bseqs + first_part + i, NULL, stdout);
    }
    printf("\nAfter appending the second part, existing graph has %d sequences, %d nodes. Now append the remaining %d sequences together.\n\n", ab->abs->n_seq, ab->abg->node_n, n_seqs - first_part - second_part);
    // 4. append the third part of sequences together
    abpt->out_msa = 1, abpt->out_cons = 1;
    abpoa_msa(ab, abpt, third_part, NULL, seq_lens + first_part + second_part, bseqs + first_part + second_part, NULL, stdout);

    printf("\nFinal graph has %d sequences, %d nodes.\n", ab->abs->n_seq, ab->abg->node_n);

    // free seq-related variables
    for (i = 0; i < n_seqs; ++i) { free(bseqs[i]); }
    free(bseqs); free(seq_lens);

    // free abpoa-related variables
    abpoa_free(ab); abpoa_free_para(abpt); 
    return 0;
}
