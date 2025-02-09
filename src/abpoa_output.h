#ifndef ABPOA_OUTPUT_H
#define ABPOA_OUTPUT_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cand_het_pos_t {
    int pos, depth, var_type;
    int n_uniq_alles;
    uint8_t *alle_bases;
} cand_het_pos_t;

typedef struct read_het_profile_t {
    int read_id; // 0 .. n_seq-1
    int start_het_idx, end_het_idx; // 0 .. n_het_pos-1
    int *alleles; // 0:ref(het->base), 1:alt
} read_het_profile_t;

typedef struct cand_het_t {
    // static information
    int pos, phase_set; // 0-based
    int var_type; // BAM_CINS/BAM_CDEL/BAM_CDIFF
    int n_depth; // may not be n_seq, as some reads are partially aligned
    int n_uniq_alles; 
    int *alle_covs; // size: n_uniq_alles
    uint8_t *alle_bases; // size: n_uniq_alles

    // dynamic information, update during haplotype assignment
    uint8_t *alle_to_hap; // var-wise (alle_to_hap): alle_i -> 1:H1/2:H2/0:not set yet
    int **hap_to_alle_profile; // read-wise: 1:H1/2:H2 -> alle_i -> read count
    int *hap_to_cons_alle; // HAP-wise (hap_to_cons_alle_i): 1:H1/2:H2 -> alle_i

    uint8_t is_low_qual, is_skipped;
} cand_het_t;


void set_65536_table(void);
void set_bit_table16(void);

#ifdef __cplusplus
}
#endif

#endif