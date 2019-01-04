/* example.c libabpoa usage example */
/* To compile: */
/*     gcc -O2 example.c -I ./include -L ./lib -labpoa -lz -o example */
/* Or: gcc -O2 example.c -I ./include ./lib/libabpoa.a -lz -o example */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "abpoa.h"

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
    int i, j, seq_n = 6;
    char seqs[6][100] = {
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT" };

    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // score parameters
    abpt->m = 5;
    abpt->match = 5;    // match score
    abpt->mismatch = 4; // mismatch penalty
    abpt->gap_open = 0; // gap open penalty
    abpt->gap_ext = 8;  // gap extension penalty 
                        // gap_penalty = gap_open + gap_len * gap_ext
    abpt->mat = (int*)malloc(abpt->m * abpt->m * sizeof(int));
    gen_simple_mat(abpt->m, abpt->mat, abpt->match, abpt->mismatch);
    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    abpt->out_pog = 1; // generate parital order graph using DOT, set 0 to disable

    int *seq_node_ids_l = (int*)malloc(seq_n * sizeof(int)); int **seq_node_ids = (int**)malloc(seq_n * sizeof(int*));

    for (i = 0; i < seq_n; ++i) {
        char *seq1 = seqs[i]; int seq_len = strlen(seq1);
        uint8_t *bseq = (uint8_t*)malloc(seq_len * sizeof(uint8_t));
        for (j = 0; j < seq_len; ++j) bseq[j] = nst_nt4_table[(int)seq1[j]];

        seq_node_ids[i] = (int*)malloc(seq_len * sizeof(int));

        abpoa_cigar_t *abpoa_cigar = 0; int n_cigar = 0;
        abpoa_align_sequence_with_graph(ab, bseq, seq_len, abpt, &n_cigar, &abpoa_cigar);

        seq_node_ids_l[i] = 0;
        abpoa_add_graph_alignment(ab->abg, abpt, bseq, seq_len, n_cigar, abpoa_cigar, seq_node_ids[i], seq_node_ids_l+i);

        free(bseq);
        if (n_cigar) free(abpoa_cigar);
    }
    /* generate consensus sequence from graph */
    if (abpt->out_cons && ab->abg->node_n > 2) {
        abpoa_generate_consensus(ab->abg, abpt->cons_agrm);
        fprintf(stdout, ">Consensus_sequence\n");
        for (i = 0; i < ab->abg->cons_l; ++i) {
            fprintf(stdout, "%c", "ACGTN"[ab->abg->cons_seq[i]]);
        } fprintf(stdout, "\n");
    }

    /* generate multiple sequence alignment */
    if (abpt->out_msa &&  ab->abg->node_n > 2)
        abpoa_generate_multiple_sequence_alingment(ab->abg, seq_node_ids, seq_node_ids_l, seq_n, abpt->out_cons, stdout);

    /* generate DOT partial order graph plot */
    if (abpt->out_pog)
        abpoa_graph_visual(ab->abg, "abpoa.dot");

    for (i = 0; i < seq_n; ++i) free(seq_node_ids[i]); free(seq_node_ids); free(seq_node_ids_l);
    abpoa_free_para(abpt); abpoa_free(ab);
    return 0;
}
