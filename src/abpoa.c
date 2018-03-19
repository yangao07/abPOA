#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "abpoa.h"
#include "abpoa_graph.h"
#include "abpoa_align.h"
#include "abpoa_graph_visual.h"
#include "align.h"
#include "simd_instruction.h"
#include "utils.h"

abpoa_para_t *abpoa_init_para(void) {
    abpoa_para_t *abpt = (abpoa_para_t*)_err_malloc(sizeof(abpoa_para_t));
    abpt->align_mode = POA_GLOBAL_FLAG;
    abpt->use_ada = 1; // use adaptive band

    // number of residue types
    abpt->m = 5; // nucleotides
    abpt->mat = NULL; // TODO score matrix for aa

    // score matrix
    abpt->match = POA_MATCH;
    abpt->mismatch = POA_MISMATCH;
    abpt->gap_open = POA_GAP_OPEN;
    abpt->gap_ext = POA_GAP_EXT;

    abpt->inf_min = MAX_OF_TWO(INF_32_MIN + POA_MISMATCH, INF_32_MIN + POA_GAP_OPEN + POA_GAP_EXT);

    abpt->bw = 10; // TODO band width
    abpt->zdrop = 100;
    abpt->end_bonus = 5;

    abpt->simd_flag = simd_check();

    return abpt;
}

void abpoa_free_para(abpoa_para_t *abpt) {
    if (abpt->mat != NULL) free(abpt->mat);
    free(abpt);
}

void cons_to_seq_score(int cons_l, uint8_t *cons_seq, int seq_n, char (*seq)[100], abpoa_para_t *abpt) {
    int i, j;
    abpt->bw = -1; // disable band width
    printf("cons_to_seq:\n");
    uint8_t *bseq = (uint8_t*)_err_malloc(100 * sizeof(uint8_t)); int bseq_m = 100;
    // progressively partial order alignment
    for (i = 0; i < seq_n; ++i) {
        abpoa_graph_t *graph = abpoa_init_graph(1);
        int seq_l = strlen(seq[i]);
        if (seq_l > bseq_m) {
            bseq_m = seq_l;
            bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
        }
        abpoa_cigar_t *abpoa_cigar=0; int n_cigar=0;
        abpoa_align_sequence_with_graph(graph, cons_seq, cons_l, abpt, &n_cigar, &abpoa_cigar); 
        abpoa_add_graph_alignment(graph, cons_seq, cons_l, n_cigar, abpoa_cigar, NULL, NULL); if (n_cigar) free(abpoa_cigar);

        for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
        abpoa_align_sequence_with_graph(graph, bseq, seq_l, abpt, &n_cigar, &abpoa_cigar); if (n_cigar) free(abpoa_cigar);
        abpoa_free_graph(graph, 1);
    }

    int rest_n = 0, rest_l, rest_i;
    char rest[100][100] = { 
        "GTCGTAAAGAACGTAGGTCGCCCGTCCGTAATCTGTCGGATCACCGGAAAGATGACGACCCGTAAAGTGATAATGATCAT",
        "CATAAAGAACGTAGGTCGCCGTGAGTCCGTAATCCGTACGGATTCACCGGAATGGCGTAGTTACCCGATAAAGTGATAATACAT",
    };
    uint8_t rest_seq[100];
    for (rest_i = 0; rest_i < rest_n; ++rest_i) {
        printf("REST# %d: cons_to_seq:\n", rest_i);
        rest_l = strlen(rest[rest_i]);
        for (i = 0; i < rest_l; ++i) rest_seq[i] = nst_nt4_table[(int)(rest[rest_i][i])];
        // progressively partial order alignment
        for (i = 0; i < seq_n; ++i) {
            abpoa_graph_t *graph = abpoa_init_graph(1);
            int seq_l = strlen(seq[i]);
            if (seq_l > bseq_m) {
                bseq_m = seq_l;
                bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
            }
            abpoa_cigar_t *abpoa_cigar=0; int n_cigar=0;
            abpoa_align_sequence_with_graph(graph, rest_seq, rest_l, abpt, &n_cigar, &abpoa_cigar); 
            abpoa_add_graph_alignment(graph, rest_seq, rest_l, n_cigar, abpoa_cigar, NULL, NULL); if (n_cigar) free(abpoa_cigar);

            for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
            abpoa_align_sequence_with_graph(graph, bseq, seq_l, abpt, &n_cigar, &abpoa_cigar); if (n_cigar) free(abpoa_cigar);
            abpoa_free_graph(graph, 1);
        }
    }
    free(bseq);
}

// int abpoa_main(const char *seq_fn, abpoa_para_t *abpt) { TODO
int abpoa_main(int seq_n, char (*seq)[100], abpoa_para_t *abpt){
    int i, j;
    abpoa_graph_t *graph = abpoa_init_graph(1);

    int **seq_node_ids = (int**)_err_malloc(seq_n * sizeof(int*));
    int *seq_node_ids_l = (int*)_err_malloc(seq_n * sizeof(int));
    uint8_t *bseq = (uint8_t*)_err_malloc(100 * sizeof(uint8_t)); int bseq_m = 100;
    // progressively partial order alignment
    for (i = 0; i < seq_n; ++i) {
        int seq_l = strlen(seq[i]);
#ifdef __DEBUG__
        printf("seq(%d): %s\n", seq_l, seq[i]);
#endif
        if (seq_l > bseq_m) {
            bseq_m = seq_l;
            bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
        }
        seq_node_ids[i] = (int*)_err_malloc(seq_l * sizeof(int));
        for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
        abpoa_cigar_t *abpoa_cigar=0; int n_cigar=0;
        abpoa_align_sequence_with_graph(graph, bseq, seq_l, abpt, &n_cigar, &abpoa_cigar);
        seq_node_ids_l[i] = 0;
        abpoa_add_graph_alignment(graph, bseq, seq_l, n_cigar, abpoa_cigar, seq_node_ids[i], seq_node_ids_l+i);
        //printf("%d %d %d\n", seq_l, graph->rank_n, seq_node_ids_l[i]);
        char abpoa_dot_fn[100]; sprintf(abpoa_dot_fn, "./dot_plot/abpoa_%d.dot", i);
        //abpoa_graph_visual(graph, abpoa_dot_fn);
        if (n_cigar) free(abpoa_cigar);
    }
    for (i = 0; i < seq_n; ++i) {
        printf("%s\n", seq[i]);
    }

    // generate consensus from graph
    abpoa_generate_consensus(graph);
    /*printf("consensus:\n");
    for (i = 0; i < graph->cons_l; ++i) {
        printf("%c", "ACGTN"[graph->cons_seq[i]]);
    } printf("\n");*/
    // generate multiple sequence alignment
    abpoa_generate_multiple_sequence_alingment(graph, seq_node_ids, seq_node_ids_l, seq_n, 1, stdout);
    //cons_to_seq_score(graph->cons_l, graph->cons_seq, seq_n, seq, abpt);

    abpoa_free_graph(graph, 1); free(bseq);
    for (i = 0; i < seq_n; ++i) free(seq_node_ids[i]); free(seq_node_ids); free(seq_node_ids_l);
    return 0;
}

int main(int argc, char **argv) {
    int seq_n = 6;
    char seq[100][100] = {
        //"ACGAATAG",
        //"ACGTAG",
        //"ATCAG",
        //"ATCAG",
        //"ATCAG",
        //"AGTAG",
        //"AGTAG",
        //"AGTCG",
        //"AGTCG",
        //"AGTCG",
        //"AGTCG",
        //"AAAATCGGCCCC",
        //"AAAATCGGCCCC",
        //"AAAATCGGCCCC",
        //"CCCCAGCATTTT",
        //"CCCCAGCATTTT",
        //"GGGGCTACTTTT",
        //"GGGGCTACTTTT",
        //"AACAT",
        //"AACT",
        //"AACT"
        //"AACAT",
        //"ACT",
        //"TTACGGAAGG",
        //"ACGGAA",
        //"ACGTAA",
        //"ACGGTTAA",
        //"ACTTGGAA",
        //"AGGGTTAA",
        //"AGTACCAA",
        //"AGGGCCAA",
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT",
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
    };
    abpoa_para_t *abpt = abpoa_init_para();
    // parse argv TODO

    abpt->mat = (int*)_err_malloc(abpt->m * abpt->m * sizeof(int));
    gen_simple_mat(abpt->m, abpt->mat, abpt->match, abpt->mismatch);
    abpoa_main(seq_n, seq, abpt);
    abpoa_free_para(abpt);

    return 0;
}
