/* To compile: 
   gcc -O3 msa_abPOA.c -I ../include -L ../lib -labpoa -lz -o msa_abPOA 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <getopt.h>
#include "abpoa.h"

#define MAX_SEQ_LEN 20000
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

char *get_seq_for_racon(FILE *fp, int *start) {
    char c = fgetc(fp);
    if (feof(fp)) return NULL; // end of file
    else if (c != '>') {
        fprintf(stderr, "Unexpected format.\n");
        exit(1);
    }
    // starts with '>'
    char *seq = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    fgets(seq, MAX_SEQ_LEN, fp); // header line
    if (seq[0] == '0') *start = 1;
    else *start = 0;
    if(fgets(seq, MAX_SEQ_LEN, fp) == NULL) { // sequence line
        fprintf(stderr, "No fasta sequence.\n");
        exit(1); 
    }

    int seq_len = strlen(seq) - 1;
    seq[seq_len] = '\0';
    return seq;
}

// read fasta format
char *get_seq(FILE *fp)
{
    char c = fgetc(fp);
    if (feof(fp)) return NULL; // end of file
    else if (c != '>') {
        fprintf(stderr, "Unexpected format.\n");
        exit(1);
    }
    // starts with '>'
    char *seq = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    fgets(seq, MAX_SEQ_LEN, fp); // header line
    if(fgets(seq, MAX_SEQ_LEN, fp) == NULL) { // sequence line
        fprintf(stderr, "No fasta sequence.\n");
        exit(1); 
    }

    int seq_len = strlen(seq) - 1;
    seq[seq_len] = '\0';
    return seq;
}

void help() {
    printf(
        "\n"
        "usage: ./msa_abPOA [options ...]\n"
        "\n"
        "    options:\n"
        "        -b <int>\n"
        "            default: 10\n"
        "            1st part of extra-band-width, set as negative to disable adaptive banded DP\n"
        "        -m <int>\n"
        "            default: 2\n"
        "            score for matching bases\n"
        "        -x <int>\n"
        "            default: 4\n"
        "            penalty for mismatching bases\n"
        "        -o <int(,int)>\n"
        "            default: gap_open1 = 4, gap_open2 = 24\n"
        "            gap opening penalty (must be non-negative)\n"
        "        -e <int(,int)>\n"
        "            default: gap_ext1 = 2, gap_ext2 = 1\n"
        "            gap extension penalty (must be non-negative)\n"
        "        -s <file>\n"
        "            default: seq.fa\n"
        "            the input sequence set\n"
        "        -n <int>\n"
        "            default: 10\n"
        "            number of sequences in each set\n"
        "        -r \n"
        "            input sequence set is from racon\n"
        "        -h \n"
        "            prints the usage\n"
    );
}

int main (int argc, char * const argv[]) {
    char *seq_file = "seq.fa";
    //initial
    abpoa_para_t *abpt = abpoa_init_para();

    int n_seqs = 0, *seq_lens; uint8_t **bseqs;
    char opt; char *s;
    int for_racon = 0;
    while ((opt = getopt(argc, argv, "l:m:x:o:e:s:n:rhb:")) != -1) {
        switch (opt) {
            case 'l': abpt->align_mode = atoi(optarg); break;
            case 'm': abpt->match = atoi(optarg); break;
            case 'x': abpt->mismatch = atoi(optarg); break;
            case 'o': abpt->gap_open1 = strtol(optarg, &s, 10); if (*s == ',') abpt->gap_open2 = strtol(s+1, &s, 10); break;
            case 'e': abpt->gap_ext1 = strtol(optarg, &s, 10); if (*s == ',') abpt->gap_ext2 = strtol(s+1, &s, 10); break;
            case 's': seq_file = optarg; break;
            case 'n': n_seqs = atoi(optarg); break;
            case 'r': for_racon = 1; break;
            case 'b': abpt->wb = atoi(optarg); break;
            case 'h': help(); return 0;
            default: help(); return 1;
        }
    }


    // output options
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable

    abpoa_post_set_para(abpt);

    FILE *fp_seq;
    if((fp_seq = fopen(seq_file, "r")) == NULL) {
           printf("fail to read seq\n");
           exit(1);
    }

    struct timeval start_time, end_time;
    double runtime = 0; int i; char *seq;
    abpoa_t *ab = abpoa_init();
    if (for_racon) {
        n_seqs = 500;
        seq_lens = (int*)malloc(sizeof(int) * n_seqs);
        bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
        int start, init=1, seq_i = 0;
        while (1) {
            while (1) {
                seq = get_seq_for_racon(fp_seq, &start);
                if (seq == NULL || start == 1) break;
                int seq_len = strlen(seq);
                seq_lens[seq_i] = seq_len;
                bseqs[seq_i] = (uint8_t*)malloc(seq_len * sizeof(uint8_t));
                for (i = 0; i < seq_len; ++i)
                    bseqs[seq_i][i] = nst_nt4_table[(int)seq[i]];
                free(seq); 
                seq_i++;
            }
            if (init) init = 0;
            else {
                // do msa
                // printf("%d\n", seq_i);
                gettimeofday(&start_time, NULL);
                abpoa_msa(ab, abpt, seq_i, seq_lens, bseqs, stdout, NULL, NULL, NULL, NULL, NULL);
                gettimeofday(&end_time, NULL);
                runtime = runtime + (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) * 1e-6;

                if (seq_i > 0) {
                    for (i = 0; i < seq_i; ++i) free(bseqs[i]); 
                }
                if(feof(fp_seq)) break;
            }
            if (start) {
                int seq_len = strlen(seq);
                seq_lens[0] = seq_len;
                bseqs[0] = (uint8_t*)malloc(seq_len * sizeof(uint8_t));
                for (i = 0; i < seq_len; ++i)
                    bseqs[0][i] = nst_nt4_table[(int)seq[i]];
                free(seq); 
                seq_i = 1;
            }
        }
    } else {
        seq_lens = (int*)malloc(sizeof(int) * n_seqs);
        bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
        while(1) {
            int seq_i = 0;
            while(1) {
                seq = get_seq(fp_seq);
                if(seq == NULL ) break;

                int seq_len = strlen(seq);
                seq_lens[seq_i] = seq_len;
                bseqs[seq_i] = (uint8_t*)malloc(seq_len * sizeof(uint8_t));
                for (i = 0; i < seq_len; ++i)
                    bseqs[seq_i][i] = nst_nt4_table[(int)seq[i]];
                free(seq); 
                seq_i++;
                if (seq_i == n_seqs) break; // read a set of sequences
            }
            if (seq_i == 0) break;

            gettimeofday(&start_time, NULL);
            abpoa_msa(ab, abpt, n_seqs, seq_lens, bseqs, stdout, NULL, NULL, NULL, NULL, NULL);
            gettimeofday(&end_time, NULL);
            runtime = runtime + (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) * 1e-6;

            if (n_seqs > 0) {
                for (i = 0; i < n_seqs; ++i) free(bseqs[i]); 
            }

            if(feof(fp_seq)) break;
        }
    }
    abpoa_free(ab, abpt);
    fprintf(stderr, "%.2f ", runtime);

    free(seq_lens); free(bseqs);
    abpoa_free_para(abpt); fclose(fp_seq);
    return 0;
}
