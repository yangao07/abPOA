//TODO: check these are needed
#include <stdio.h>

#include "defines.h"
#include "structures.h"

/////////////////////////////////////////////////////////////
// Function: print_sim_param - Print simulation parameters //
/////////////////////////////////////////////////////////////

void print_sim_param(const sim_t *sim) {
  fprintf(stderr, ":::: Simulation parameters :::\n\n");

  if (sim->process == PROCESS_MODEL) {
    fprintf(stderr, "Simulated by stochastic model.\n\n");
  } else {
    fprintf(stderr, "Simulated by fastq sampling.\n\n");
  }

  fprintf(stderr, "prefix : %s\n", sim->prefix);
  if (sim->set_flg[14]) {
    fprintf(stderr, "sample_profile_id : %s\n", sim->profile_id);
  }

  if (sim->data_type == DATA_TYPE_CLR) {
    fprintf(stderr, "data-type : CLR\n");
  } else {
    fprintf(stderr, "data-type : CCS\n");
  }

  fprintf(stderr, "depth : %lf\n", sim->depth);

  if (sim->set_flg[0]) {
    fprintf(stderr, "length-mean : (sampling FASTQ)\n");
    fprintf(stderr, "length-sd : (sampling FASTQ)\n");
  } else {
    fprintf(stderr, "length-mean : %f\n", sim->len_mean);
    fprintf(stderr, "length-sd : %f\n", sim->len_sd);
  }
  fprintf(stderr, "length-min : %ld\n", sim->len_min);
  fprintf(stderr, "length-max : %ld\n", sim->len_max);

  if (sim->set_flg[0]) {
    fprintf(stderr, "accuracy-mean : (sampling FASTQ)\n");
    fprintf(stderr, "accuracy-sd : (sampling FASTQ)\n");
  } else {
    fprintf(stderr, "accuracy-mean : %f\n", sim->accuracy_mean);
    fprintf(stderr, "accuracy-sd : %f\n", sim->accuracy_sd);
  }
  fprintf(stderr, "accuracy-min : %f\n", sim->accuracy_min);
  fprintf(stderr, "accuracy-max : %f\n", sim->accuracy_max);

  fprintf(stderr, "difference-ratio : %ld:%ld:%ld\n",
    sim->sub_ratio, sim->ins_ratio, sim->del_ratio);

  fprintf(stderr, "\n");
}

/////////////////////////////////////////////////////
// Function: print_fastq_stats - Print FASTQ stats //
/////////////////////////////////////////////////////

void print_fastq_stats(const sim_t *sim, const fastq_t *fastq) {
  fprintf(stderr, ":::: FASTQ stats ::::\n\n");

  if (sim->process == PROCESS_SAMPLING_REUSE) {
    fprintf(stderr, "file name : %s\n", sim->profile_fq);
  } else {
    fprintf(stderr, "file name : %s\n", fastq->file);
    fprintf(stderr, "\n:: all reads ::\n");
    fprintf(stderr, "read num. : %ld\n", fastq->num);
    fprintf(stderr, "read total length : %lld\n", fastq->len_total);
    fprintf(stderr, "read min length : %ld\n", fastq->len_min);
    fprintf(stderr, "read max length : %ld\n", fastq->len_max);
  }

  fprintf(stderr, "\n:: filtered reads ::\n");
  fprintf(stderr, "read num. : %ld\n", fastq->num_filtered);
  fprintf(stderr, "read total length : %lld\n", fastq->len_total_filtered);
  fprintf(stderr, "read min length : %ld\n", fastq->len_min_filtered);
  fprintf(stderr, "read max length : %ld\n", fastq->len_max_filtered);
  fprintf(stderr, "read length mean (SD) : %f (%f)\n",
    fastq->len_mean_filtered, fastq->len_sd_filtered);
  fprintf(stderr, "read accuracy mean (SD) : %f (%f)\n",
    fastq->accuracy_mean_filtered, fastq->accuracy_sd_filtered);
  fprintf(stderr, "\n");
}

////////////////////////////////////////////////////////////////
// Function: print_simulation_stats - Print Simulation Stats. //
////////////////////////////////////////////////////////////////

void print_simulation_stats(const sim_t *sim, const ref_t *ref) {
  double res_depth = (double)sim->res_len_total / ref->len;
  double res_sub_rate = (double)sim->res_sub_num / sim->res_len_total;
  double res_ins_rate = (double)sim->res_ins_num / sim->res_len_total;
  double res_del_rate = (double)sim->res_del_num / sim->res_len_total;

  fprintf(stderr, ":::: Simulation stats (ref.%ld - name: %s) ::::\n\n", ref->num, ref->id);
  fprintf(stderr, "read num. : %ld\n", sim->res_num);
  fprintf(stderr, "depth : %lf\n", res_depth);
  fprintf(stderr, "read length mean (SD) : %f (%f)\n",
    sim->res_len_mean, sim->res_len_sd);
  fprintf(stderr, "read length min : %ld\n", sim->res_len_min);
  fprintf(stderr, "read length max : %ld\n", sim->res_len_max);
  fprintf(stderr, "read accuracy mean (SD) : %f (%f)\n",
    sim->res_accuracy_mean, sim->res_accuracy_sd);
  fprintf(stderr, "substitution rate. : %f\n", res_sub_rate);
  fprintf(stderr, "insertion rate. : %f\n", res_ins_rate);
  fprintf(stderr, "deletion rate. : %f\n", res_del_rate);
  fprintf(stderr, "\n");
}

///////////////////////////////////////
// Function: print_help - Print help //
///////////////////////////////////////

void print_help() {
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: pbsim [options] <reference>\n\n");
  fprintf(stderr, " <reference>           FASTA format file.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [general options]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --prefix             prefix of output files (sd).\n");
  fprintf(stderr, "  --data-type          data type. CLR or CCS (CLR).\n");
  fprintf(stderr, "  --depth              depth of coverage (CLR: 20.0, CCS: 50.0).\n");
  fprintf(stderr, "  --length-min         minimum length (100).\n");
  fprintf(stderr, "  --length-max         maximum length (CLR: 25000, CCS: 2500).\n");
  fprintf(stderr, "  --accuracy-min       minimum accuracy.\n");
  fprintf(stderr, "                       (CLR: 0.75, CCS: fixed as 0.75).\n");
  fprintf(stderr, "                       this option can be used only in case of CLR.\n");
  fprintf(stderr, "  --accuracy-max       maximum accuracy.\n");
  fprintf(stderr, "                       (CLR: 1.00, CCS: fixed as 1.00).\n");
  fprintf(stderr, "                       this option can be used only in case of CLR.\n");
  fprintf(stderr, "  --difference-ratio   ratio of differences. substitution:insertion:deletion.\n");
  fprintf(stderr, "                       each value up to 1000 (CLR: 10:60:30, CCS:6:21:73).\n");
  fprintf(stderr, "  --seed               for a pseudorandom number generator (Unix time).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options of sampling-based simulation]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --sample-fastq       FASTQ format file to sample.\n");
  fprintf(stderr, "  --sample-profile-id  sample-fastq (filtered) profile ID.\n");
  fprintf(stderr, "                       when using --sample-fastq, profile is stored.\n");
  fprintf(stderr, "                       'sample_profile_<ID>.fastq', and\n");
  fprintf(stderr, "                       'sample_profile_<ID>.stats' are created.\n");
  fprintf(stderr, "                       when not using --sample-fastq, profile is re-used.\n");
  fprintf(stderr, "                       Note that when profile is used, --length-min,max,\n");
  fprintf(stderr, "                       --accuracy-min,max would be the same as the profile.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options of model-based simulation].\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --model_qc           model of quality code.\n");
  fprintf(stderr, "  --length-mean        mean of length model (CLR: 3000.0, CCS:450.0).\n");
  fprintf(stderr, "  --length-sd          standard deviation of length model.\n");
  fprintf(stderr, "                       (CLR: 2300.0, CCS: 170.0).\n");
  fprintf(stderr, "  --accuracy-mean      mean of accuracy model.\n");
  fprintf(stderr, "                       (CLR: 0.78, CCS: fixed as 0.98).\n");
  fprintf(stderr, "                       this option can be used only in case of CLR.\n");
  fprintf(stderr, "  --accuracy-sd        standard deviation of accuracy model.\n");
  fprintf(stderr, "                       (CLR: 0.02, CCS: fixed as 0.02).\n");
  fprintf(stderr, "                       this option can be used only in case of CLR.\n");
  fprintf(stderr, "\n");
}