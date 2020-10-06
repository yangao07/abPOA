/* sub_example.c libabpoa usage example
   To compile:
      gcc -O2 sub_example.c -I ./include -L ./lib -labpoa -lz -o sub_example
  Or: gcc -O2 sub_example.c -I ./include ./lib/libabpoa.a -lz -o sub_example
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "src/abpoa.h"

// AaCcGgTtNn ==> 0,1,2,3,4
unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int main(void) {
    int i, j, n_seqs = 4;
    char seqs[100][1000] = {
         // 0       1         2         3
         // 23456789012345678901234567890123               
           "CGTCAATCTATCGAAGCATACGCGGGCAGAGC",
        "CCACGTCAATCTATCGAAGCATACGCGGCAGC",
               "AATCTATCGAAGCATACG",
              "CAATGCTAGTCGAAGCAGCTGCGGCAG",
           "CGTCAATCTATCGAAGCATTCTACGCGGCAGAGC",
        "CGTCAATCTAGAAGCATACGCGGCAAGAGC",
        "CGTCAATCTATCGGTAAAGCATACGCTCTGTAGC",
        "CGTCAATCTATCTTCAAGCATACGCGGCAGAGC",
        "CGTCAATGGATCGAGTACGCGGCAGAGC",
        "CGTCAATCTAATCGAAGCATACGCGGCAGAGC"
        };

    int beg_end_id[100][2] = {
        {0, 0}, 
        {2, 33},
        {6, 23}, 
        {5, 30}, 
        {0, 1}, 
        {0, 1}, 
        {0, 1}, 
        {0, 1}, 
        {0, 1}, 
        {0, 1}, 
        //{2, 52}, 
        //{2, 52}, 
        //{2, 52}, 
        //{2, 52}, 
        //{2, 52}, 
        //{2, 52}, 
        //{2, 52}, 
        //{2, 52}, 
        //{2, 52} 
    };

    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // alignment parameters
    // abpt->align_mode = 0; // 0:global alignment, 1:extension
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

    abpoa_post_set_para(abpt);

    // collect sequence length, trasform ACGT to 0123
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = strlen(seqs[i]);
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j)
            bseqs[i][j] = nst_nt4_table[(int)seqs[i][j]];
    }

    // perform abpoa-msa
    abpoa_res_t res;
    for (i = 0; i < n_seqs; ++i) {
        res.graph_cigar = 0, res.n_cigar = 0; res.is_rc = 0;
        abpoa_align_sequence_to_subgraph(ab, abpt, beg_end_id[i][0], beg_end_id[i][1], bseqs[i], seq_lens[i], &res);
        int exc_beg, exc_end;
        if (i != 0) abpoa_subgraph_nodes(ab, beg_end_id[i][0], beg_end_id[i][1], &exc_beg, &exc_end);
        else exc_beg = 0, exc_end = 0;
        abpoa_add_subgraph_alignment(ab, abpt, exc_beg, exc_end, bseqs[i], seq_lens[i], res, i, n_seqs);
        if (res.n_cigar) free(res.graph_cigar);
    }
    if (abpt->out_cons) {
        abpoa_generate_consensus(ab, abpt, n_seqs, stdout, NULL, NULL, NULL, NULL);
        if (ab->abg->is_called_cons == 0)
            fprintf(stderr, "Warning: no consensus sequence generated.\n");
    }
    if (abpt->out_msa) {
        abpoa_generate_rc_msa(ab, abpt, NULL, NULL, n_seqs, stdout, NULL, NULL);
    }
    // abpoa_reset_graph(ab, abpt, seq_lens[0]); // reset graph before re-use

    /* generate DOT partial order graph plot */
    abpt->out_pog = strdup("example.pdf"); // dump parital order graph to file
    if (abpt->out_pog != NULL) abpoa_dump_pog(ab, abpt);
    for (i = 0; i < n_seqs; ++i) free(bseqs[i]); free(bseqs); free(seq_lens);
    abpoa_free(ab, abpt); abpoa_free_para(abpt);
    return 0;
}
