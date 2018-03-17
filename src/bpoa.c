#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bpoa.h"
#include "bpoa_graph.h"
#include "bpoa_align.h"
#include "bpoa_graph_visual.h"
#include "align.h"
#include "simd_instruction.h"
#include "utils.h"

bpoa_para_t *bpoa_init_para(void) {
    bpoa_para_t *bpt = (bpoa_para_t*)_err_malloc(sizeof(bpoa_para_t));
    bpt->align_mode = POA_GLOBAL_FLAG;

    // number of residue types
    bpt->m = 5; // nucleotides
    bpt->mat = NULL; // TODO score matrix for aa

    // score matrix
    bpt->match = POA_MATCH;
    bpt->mismatch = POA_MISMATCH;
    bpt->gap_open = POA_GAP_OPEN;
    bpt->gap_ext = POA_GAP_EXT;

    bpt->inf_min = MAX_OF_TWO(INF_32_MIN + POA_MISMATCH, INF_32_MIN + POA_GAP_OPEN + POA_GAP_EXT);

    bpt->bw = 10; // TODO band width
    bpt->zdrop = 100;
    bpt->end_bonus = 5;

    bpt->simd_flag = simd_check();

    return bpt;
}

void bpoa_free_para(bpoa_para_t *bpt) {
    if (bpt->mat != NULL) free(bpt->mat);
    free(bpt);
}

void cons_to_seq_score(int cons_l, uint8_t *cons_seq, int seq_n, char (*seq)[100], bpoa_para_t *bpt) {
    int i, j;
    bpt->bw = -1; // disable band width
    printf("cons_to_seq:\n");
    uint8_t *bseq = (uint8_t*)_err_malloc(100 * sizeof(uint8_t)); int bseq_m = 100;
    // progressively partial order alignment
    for (i = 0; i < seq_n; ++i) {
        bpoa_graph_t *graph = bpoa_init_graph(1);
        int seq_l = strlen(seq[i]);
        if (seq_l > bseq_m) {
            bseq_m = seq_l;
            bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
        }
        bpoa_cigar_t *bpoa_cigar=0; int n_cigar=0;
        bpoa_align_sequence_with_graph(graph, cons_seq, cons_l, bpt, &n_cigar, &bpoa_cigar); 
        bpoa_add_graph_alignment(graph, cons_seq, cons_l, n_cigar, bpoa_cigar, NULL, NULL); if (n_cigar) free(bpoa_cigar);

        for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
        bpoa_align_sequence_with_graph(graph, bseq, seq_l, bpt, &n_cigar, &bpoa_cigar); if (n_cigar) free(bpoa_cigar);
        bpoa_free_graph(graph, 1);
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
            bpoa_graph_t *graph = bpoa_init_graph(1);
            int seq_l = strlen(seq[i]);
            if (seq_l > bseq_m) {
                bseq_m = seq_l;
                bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
            }
            bpoa_cigar_t *bpoa_cigar=0; int n_cigar=0;
            bpoa_align_sequence_with_graph(graph, rest_seq, rest_l, bpt, &n_cigar, &bpoa_cigar); 
            bpoa_add_graph_alignment(graph, rest_seq, rest_l, n_cigar, bpoa_cigar, NULL, NULL); if (n_cigar) free(bpoa_cigar);

            for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
            bpoa_align_sequence_with_graph(graph, bseq, seq_l, bpt, &n_cigar, &bpoa_cigar); if (n_cigar) free(bpoa_cigar);
            bpoa_free_graph(graph, 1);
        }
    }
    free(bseq);
}

// int bpoa_main(const char *seq_fn, bpoa_para_t *bpt) { TODO
int bpoa_main(int seq_n, char (*seq)[100], bpoa_para_t *bpt){
    int i, j;
    bpoa_graph_t *graph = bpoa_init_graph(1);

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
        bpoa_cigar_t *bpoa_cigar=0; int n_cigar=0;
        bpoa_align_sequence_with_graph(graph, bseq, seq_l, bpt, &n_cigar, &bpoa_cigar);
        seq_node_ids_l[i] = 0;
        bpoa_add_graph_alignment(graph, bseq, seq_l, n_cigar, bpoa_cigar, seq_node_ids[i], seq_node_ids_l+i);
        //printf("%d %d %d\n", seq_l, graph->rank_n, seq_node_ids_l[i]);
        char bpoa_dot_fn[100]; sprintf(bpoa_dot_fn, "./dot_plot/bpoa_%d.dot", i);
        //bpoa_graph_visual(graph, bpoa_dot_fn);
        if (n_cigar) free(bpoa_cigar);
    }
    for (i = 0; i < seq_n; ++i) {
        printf("%s\n", seq[i]);
    }

    // generate consensus from graph
    bpoa_generate_consensus(graph);
    /*printf("consensus:\n");
    for (i = 0; i < graph->cons_l; ++i) {
        printf("%c", "ACGTN"[graph->cons_seq[i]]);
    } printf("\n");*/
    // generate multiple sequence alignment
    bpoa_generate_multiple_sequence_alingment(graph, seq_node_ids, seq_node_ids_l, seq_n, 1, stdout);
    //cons_to_seq_score(graph->cons_l, graph->cons_seq, seq_n, seq, bpt);

    bpoa_free_graph(graph, 1); free(bseq);
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
    bpoa_para_t *bpt = bpoa_init_para();
    // parse argv TODO

    bpt->mat = (int*)_err_malloc(bpt->m * bpt->m * sizeof(int));
    gen_simple_mat(bpt->m, bpt->mat, bpt->match, bpt->mismatch);
    bpoa_main(seq_n, seq, bpt);
    bpoa_free_para(bpt);

    return 0;
}
