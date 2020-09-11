#pragma once 
#include "defines.h"

// simulation parameters
struct sim_t {
  int data_type;
  int process;
  double depth;
  double accuracy_mean, accuracy_sd, accuracy_max, accuracy_min;
  long len_min, len_max; 
  double len_mean, len_sd; 
  long long len_quota;
  long sub_ratio, ins_ratio, del_ratio;
  double sub_rate, ins_rate, del_rate;
  int set_flg[20];
  long res_num;
  long long res_len_total; 
  double res_accuracy_mean, res_accuracy_sd;
  long res_len_min, res_len_max; 
  double res_len_mean, res_len_sd; 
  long res_sub_num, res_ins_num, res_del_num;
  double res_sub_rate, res_ins_rate, res_del_rate;
  char *prefix, *outfile_ref, *outfile_fq, *outfile_maf;
  char *model_qc_file;
  char *profile_id, *profile_fq, *profile_stats;
};

// FASTQ
struct fastq_t {
  char *file;
  long num;
  long long len_total;
  long len_min, len_max;
  long num_filtered;
  long long len_total_filtered;
  long len_min_filtered, len_max_filtered;
  double len_mean_filtered, len_sd_filtered;
  double accuracy_mean_filtered, accuracy_sd_filtered;
};

// Reference
struct ref_t {
  char *file;
  char *seq;
  char id[REF_ID_LEN_MAX + 1];
  long len;
  long num_seq;
  long num;
};

// Mutation
struct mut_t {
  long sub_thre[94], ins_thre[94], del_thre;
  char *sub_nt_a, *sub_nt_t, *sub_nt_g, *sub_nt_c, *sub_nt_n, *ins_nt;
  char *qc, *new_qc, *tmp_qc, *seq, *new_seq, *maf_seq, *maf_ref_seq;
  long tmp_len_max;
  char seq_strand;
  long seq_left, seq_right;
};

// Quality code
struct qc_t {
  char character;
  double prob;
};

// Quality code of model
struct model_qc_t {
  int min;
  int max;
  double prob[94];
};
