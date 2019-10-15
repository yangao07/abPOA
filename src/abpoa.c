#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "abpoa.h"
#include "abpoa_graph.h"
#include "abpoa_align.h"
#include "align.h"
#include "kseq.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

char PROG[20] = "abPOA";

const struct option abpoa_long_opt [] = {
    { "in-list", 0, NULL, 'l' },
    { "align-mode", 1, NULL, 'a' },

    { "gap-mode", 1, NULL, 'g' },
    { "match", 1, NULL, 'M' },
    { "mismatch", 1, NULL, 'x' },
    { "gap-open", 1, NULL, 'o' },
    { "gap-ext", 1, NULL, 'e' },

    { "band-width", 1, NULL, 'w' },
    { "zdrop", 1, NULL, 'z' },
    { "end-bouns", 1, NULL, 'b' },

    { "out-cons", 0, NULL, 'c' },
    { "cons-agrm", 1, NULL, 'C' },
    { "multiploid", 1, NULL, 'm', },
    { "min-frequency", 1, NULL, 'f', },
    { "out-msa", 0, NULL, 's' },

    { "out-pog", 0, NULL, 'v' },

    { 0, 0, 0, 0}
};

int abpoa_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s [option] <in.fa/fq> > msa.out\n\n", PROG);
    err_printf("Options:\n\n");
    err_printf("         -l --in-list                use list of filename as input. [False]\n\n");

    err_printf("         -a --aln-mode   [INT]       align mode. [%d]\n", ABPOA_GLOBAL_MODE);
    err_printf("                                       %d: global\n", ABPOA_GLOBAL_MODE);
    err_printf("                                       %d: extension\n\n", ABPOA_EXTEND_MODE);
    //err_printf("                                       %d: local\n", ABPOA_LOCAL_MODE);
    //err_printf("                                       %d: semi-global\n\n", ABPOA_SEMI_MODE);
    err_printf("         -g --gap-mode   [INT]       gap penalty mode. [%d]\n", ABPOA_CONVEX_GAP);
    err_printf("                                     When \'-o\' and \'-e\' are set, gap mode will be changed accordingly.\n");
    err_printf("                                     See \'-o\' and \'-e\' for more details.\n");
    err_printf("                                       %d: linear gap, g*e1\n", ABPOA_LINEAR_GAP);
    err_printf("                                       %d: affine gap, o1+g*e1\n", ABPOA_AFFINE_GAP);
    err_printf("                                       %d: convex gap, min{o1+g*e1, o2+g*e2}\n\n", ABPOA_CONVEX_GAP);

    err_printf("         -M --match      [INT]       match score. [%d]\n", ABPOA_MATCH);
    err_printf("         -x --mismatch   [INT]       mismatch penalty. [%d]\n", ABPOA_MISMATCH);
    err_printf("         -o --gap-open   [INT(,INT)] gap open penalty. [%d,%d]\n", ABPOA_GAP_OPEN1, ABPOA_GAP_OPEN2);
    // TODO remind to use valid score para
    err_printf("         -e --gap-ext    [INT(,INT)] gap extension penalty [%d,%d]\n", ABPOA_GAP_EXT1, ABPOA_GAP_EXT2);
    err_printf("                                     1. abPOA uses 2-piece affine gap penalty by default, i.e.,\n");
    err_printf("                                        2-piece penalty of a g-long gap: min{o1+g*e1, o2+g*e2}\n");
    err_printf("                                     2. Set o2 as 0 to apply affine gap penalty and only o1,e1 will be used, i.e.,\n");
    err_printf("                                        affine penalty of a g-long gap: o1+g*e1\n");
    err_printf("                                     3. Set o1 as 0 to apply linear gap penalty and only e1 will be used, i.e.,\n");
    err_printf("                                        linear penalty of a g-long gap: g*e1\n\n");

    err_printf("         -w --band-width [INT]       band width used in alignment. [-1]\n");
    err_printf("                                     Effective for global and extension alignment. Set negative to disable.\n");
    //err_printf("         -a --ada-band               adaptively update band width during alignment. [False]\n");
    //err_printf("                                     Effective for global and extension alignment.\n");
    err_printf("         -z --zdrop      [INT]       Z-drop score. Effective for extension alignment. Set negative to disable. [-1]\n");
    err_printf("         -e --end-bonus  [INT]       end bonus score. Effective for extension alignment. Set negative to disable. [-1]\n\n");

    err_printf("         -c --out-cons               output consensus sequence. [False]\n");
    err_printf("         -C --cons-agrm  [INT]       algorithm for consensus calling. [0]\n");
    err_printf("                                       0: heaviest bundling\n");
    err_printf("                                       1: minimum flow\n");
    err_printf("                                       2: row-column MSA\n\n");
    err_printf("         -m --multiploid             maximum number of output consensus sequences (for multiploid data). [%d]\n", ABPOA_MULTIP);
    err_printf("                                     When 2 or more consensus are needed, --cons-agrm is set to \'row-column MSA\'\n");
    err_printf("         -f --min-frequency [DOU]    minimum frequency of each haploid (for multiploid data). [%.2f]\n", ABPOA_MIN_FRE);
    err_printf("         -s --out-msa                output multiple sequence alignment in pir format. [False]\n");
    err_printf("         -v --out-pog                generate visualized partial-order graph. [False]\n");

    err_printf("\n");
    return 1;
}

int abpoa_read_seq(kseq_t *read_seq, int chunk_read_n)
{
    kseq_t *s = read_seq;
    int n = 0;
    while (kseq_read(s+n) >= 0) {
        n++;
        if (n >= chunk_read_n) break;
    }
    return n;
}

#define abpoa_core(read_fn) {   \
    gzFile readfp = xzopen(read_fn, "r"); kstream_t *fs = ks_init(readfp);  \
    for (i = 0; i < CHUNK_READ_N; ++i) read_seq[i].f = fs;  \
    /* progressively partial order alignment */     \
    n_seqs = 0,  tot_n = 0, read_id = 0; \
    /* reset graph for a new input file */  \
    abpoa_reset_graph(ab, bseq_m, abpt);    \
    while ((n_seqs = abpoa_read_seq(read_seq, CHUNK_READ_N)) != 0) {    \
        tot_n += n_seqs;    \
        if (abpt->use_read_ids) read_ids_n = (tot_n-1) / 64 + 1;   \
        for (i = 0; i < n_seqs; ++i) {  \
            kseq_t *seq = read_seq + i; \
            int seq_l = seq->seq.l; \
            if (seq_l <= 0) err_fatal("read_seq", "Unexpected read length: %d (%s)", seq_l, seq->name.s);   \
            char *seq1 = seq->seq.s;    \
            /* printf("seq(%d): %s\n", seq_l, seq1); */   \
            if (seq_l > bseq_m) {   \
                bseq_m = MAX_OF_TWO(seq_l+1, bseq_m * 2);   \
                bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));  \
            }   \
            for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq1[j])];    \
            abpoa_cigar_t *abpoa_cigar=0; int n_cigar=0;    \
            abpoa_align_sequence_with_graph(ab, bseq, seq_l, abpt, &n_cigar, &abpoa_cigar); \
            abpoa_add_graph_alignment(ab->abg, abpt, bseq, seq_l, n_cigar, abpoa_cigar, read_id++, read_ids_n);     \
            if (n_cigar) free(abpoa_cigar); \
            /*if (abpt->out_pog) {    \
                char abpoa_dot_fn[100]; sprintf(abpoa_dot_fn, "./abpoa_%d.dot", i);    \
                abpoa_graph_visual(ab->abg, abpoa_dot_fn);  \
            }*/   \
        }   \
    }   \
    /* generate consensus from graph */ \
    if (abpt->out_cons && ab->abg->node_n > 2) {   \
        abpoa_generate_consensus(ab->abg, abpt->cons_agrm, abpt->multip, abpt->min_fre, tot_n, stdout, NULL, NULL, NULL); \
    }   \
    /* generate multiple sequence alignment */  \
    if (abpt->out_msa &&  ab->abg->node_n > 2)  \
        abpoa_generate_multiple_sequence_alingment(ab->abg, tot_n, stdout);   \
    /* generate dot plot */     \
    if (abpt->out_pog) {    \
        char abpoa_dot_fn[100] = "./abpoa.dot";    \
        abpoa_graph_visual(ab->abg, abpoa_dot_fn);  \
    }   \
    ks_destroy(fs); gzclose(readfp);    \
}

int abpoa_main(const char *list_fn, int in_list, abpoa_para_t *abpt){
    kseq_t *read_seq = (kseq_t*)calloc(CHUNK_READ_N, sizeof(kseq_t));
    int bseq_m = 1000; uint8_t *bseq = (uint8_t*)_err_malloc(bseq_m * sizeof(uint8_t));
    // int nseqs_m, *seq_m=NULL, *seq_node_ids_l=NULL, **seq_node_ids=NULL; // for msa
    int i, j, n_seqs, tot_n, read_id;
    int read_ids_n = 0;
    if (abpt->cons_agrm == ABPOA_RC || abpt->out_msa || abpt->multip > 1) {
        abpt->use_read_ids = 1;
        set_65536_table();
        if (abpt->cons_agrm == ABPOA_RC || abpt->multip > 1) set_bit_table16();
    }

    // TODO abpoa_init for each input file ???
    abpoa_t *ab = abpoa_init();
    if (in_list) { // input file list
        FILE *list_fp = fopen(list_fn, "r"); char read_fn[1024];
        while (fgets(read_fn, sizeof(read_fn), list_fp)) {
            read_fn[strlen(read_fn)-1] = '\0';
            abpoa_core(read_fn);
        }
        fclose(list_fp);
    } else { // input file
        abpoa_core(list_fn);
    }

    free(bseq);
    for (i = 0; i < CHUNK_READ_N; ++i) { free((read_seq+i)->name.s); free((read_seq+i)->comment.s); free((read_seq+i)->seq.s); free((read_seq+i)->qual.s); } free(read_seq); 
    abpoa_free(ab, abpt);
    return 0;
}

int main(int argc, char **argv) {
    int c, m, g, in_list=0, manual_gap_mode=0; char *s; abpoa_para_t *abpt = abpoa_init_para();
    while ((c = getopt_long(argc, argv, "la:g:w:z:b:M:x:o:e:m:f:csC:v", abpoa_long_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'l': in_list = 1; break;
            case 'a': m = atoi(optarg);
                      if (m != ABPOA_GLOBAL_MODE && m != ABPOA_EXTEND_MODE) { err_printf("Unknown alignment mode: %d.\n", m); return abpoa_usage(); }
                      abpt->align_mode=m; break;
            case 'g': g = atoi(optarg);
                      if (g != ABPOA_CONVEX_GAP && g != ABPOA_AFFINE_GAP && g != ABPOA_LINEAR_GAP) { err_printf("Unknown gap penalty mode: %d.\n", g); return abpoa_usage(); }
                      abpt->gap_mode = g; manual_gap_mode = 1; break;

            case 'w': abpt->bw = atoi(optarg); break;
            case 'z': abpt->zdrop = atoi(optarg); break;
            case 'b': abpt->end_bonus= atoi(optarg); break;

            case 'M': abpt->match = atoi(optarg); break;
            case 'x': abpt->mismatch = atoi(optarg); break;
            case 'o': abpt->gap_open1 = strtol(optarg, &s, 10); if (*s == ',') abpt->gap_open2 = strtol(s+1, &s, 10); break;
            case 'e': abpt->gap_ext1 = strtol(optarg, &s, 10); if (*s == ',') abpt->gap_ext2 = strtol(s+1, &s, 10); break;

            case 'c': abpt->out_cons = 1; break;
            case 'C': abpt->cons_agrm = atoi(optarg); break;
            case 'm': abpt->multip = atoi(optarg); break;
            case 'f': abpt->min_fre = atof(optarg); break;
            case 's': abpt->out_msa = 1; break;

            case 'v': abpt->out_pog= 1; break;
            default:
                      err_printf("Error: unknown option: %s.\n", optarg);
                      return abpoa_usage();
                      break;
        }
    }

    if (argc - optind != 1) return abpoa_usage();

    abpt->mat = (int*)_err_malloc(abpt->m * abpt->m * sizeof(int));
    if (manual_gap_mode == 0) abpoa_set_gap_mode(abpt);
    gen_simple_mat(abpt->m, abpt->mat, abpt->match, abpt->mismatch);
    abpoa_main(argv[optind], in_list, abpt);
    abpoa_free_para(abpt);
    return 0;
}
