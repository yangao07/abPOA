#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <algorithm> // min

#include "defines.h"
#include "structures.h"
#include "output.h"  // printing code is here
#include "helpers.h" // string manipulation and time

/////////////////////////////////////////
// Global variances                    //
/////////////////////////////////////////

FILE *fp_filtered, *fp_stats, *fp_ref, *fp_fq, *fp_maf;
struct fastq_t fastq;
struct mut_t mut;
struct model_qc_t model_qc[101];

long freq_len[FASTQ_LEN_MAX + 1];
long freq_accuracy[100000 + 1];
const char *versionString = "1.0.4";


/////////////////////////////////////////
// Prototypes of functions             //
/////////////////////////////////////////

void parse_options(int argc, char** argv, sim_t *sim);
void init_sim_res(sim_t *sim);
int initialize_sim_process(sim_t *sim, fastq_t *fastq, qc_t *qc);
int set_sim_param(sim_t *sim);
int get_ref_inf(const sim_t *sim, ref_t *ref);
int get_ref_seq(const sim_t *sim, ref_t *ref);
int get_fastq_inf(const sim_t *sim, const qc_t *qc);
int set_model_qc(const sim_t *sim);
int set_mut(sim_t *sim);
int simulate_by_sampling(sim_t *sim, ref_t *ref, mut_t *mut, fastq_t *fastq, qc_t *qc);
int simulate_by_model(sim_t *sim, ref_t *ref, qc_t *qc);
int mutate(sim_t *sim, ref_t *ref);

/////////////////////////////////////////
// Main                                //
/////////////////////////////////////////

int main (int argc, char** argv) {
  long len;
  long i, j;
  long rst1, rst2;
  long t1, t2;
  struct sim_t sim;
  struct qc_t qc[94];
  struct ref_t ref;

  rst1 = get_time_cpu();
  t1 = get_time();
  memset(sim.set_flg, 0, sizeof(sim.set_flg));

  ////// Check Input
  parse_options(argc, argv, &sim);

  if (argv[optind] == NULL) {
    print_help();
    exit(-1);
  }

  ///// Initialize

  // Quality code to error probability
  for (i=0; i<=93; i++) {
    qc[i].prob = pow(10, (double)i / -10);
    qc[i].character = (char)(i+33);
  }

  if(initialize_sim_process(&sim, &fastq, qc) == FAILED)
    exit(-1);

  // Reference sequence
  if ((ref.file = (char *)malloc(strlen(argv[optind]) + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    exit(-1);
  }
  strcpy(ref.file, argv[optind]);
  fprintf(stderr, "!! ref file %s\n", ref.file);

  if (get_ref_inf(&sim, &ref) == FAILED) {
    exit(-1);
  }

  // Set mutation parameters and varianeces
  if (set_mut(&sim) == FAILED) {
    exit(-1);
  }

  // Creating simulated reads
  for (ref.num=1; ref.num<=ref.num_seq; ref.num++) {
    if (get_ref_seq(&sim, &ref) == FAILED) {
      exit(-1);
    }

    init_sim_res(&sim);

    sprintf(sim.outfile_fq, "%s_%04ld.fastq", sim.prefix, ref.num);
    if ((fp_fq = fopen(sim.outfile_fq, "w")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_fq);
      return FAILED;
    }

    sprintf(sim.outfile_maf, "%s_%04ld.maf", sim.prefix, ref.num);
    if ((fp_maf = fopen(sim.outfile_maf, "w")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_maf);
      fclose(fp_fq);
      return FAILED;
    }

    sim.len_quota = (long long)(sim.depth * ref.len);

    if (sim.process == PROCESS_MODEL) {
      if (simulate_by_model(&sim, &ref, qc) == FAILED) {
        exit(-1);
      }
    } else {
      if (simulate_by_sampling(&sim, &ref, &mut, &fastq, qc) == FAILED) {
        exit(-1);
      }
    }

    print_simulation_stats(&sim, &ref);

    fclose(fp_fq);
    fclose(fp_maf);
  }

  if ((sim.process == PROCESS_SAMPLING_STORE) || (sim.process == PROCESS_SAMPLING_REUSE)) {
    fclose(fp_filtered);
    fclose(fp_stats);
  }

  rst2 = get_time_cpu();
  t2 = get_time();

  fprintf(stderr, ":::: System utilization ::::\n\n");
  fprintf(stderr, "CPU time(s) : %ld\n", rst2 - rst1);
  fprintf(stderr, "Elapsed time(s) : %ld\n", t2 - t1);

  return(0);
}

///////////////////////////////////////////////////////////////
// Function: initialize_sim_process - perform the sim process requested //
///////////////////////////////////////////////////////////////

int initialize_sim_process(sim_t *sim, fastq_t *fastq, qc_t *qc)
{
// Setting of simulation parameters
  if (set_sim_param(sim) == FAILED) {
    exit(-1);
  }
  print_sim_param(sim);

  switch(sim->process){
    case PROCESS_SAMPLING:
	    if ((fp_filtered = tmpfile()) == NULL) {
	      fprintf(stderr, "ERROR: Cannot open temporary file\n");
	      return FAILED;
	    }
      	break;
    case PROCESS_SAMPLING_STORE:
	     if ((fp_filtered = fopen(sim->profile_fq, "w+")) == NULL) {
	      fprintf(stderr, "ERROR: Cannot open filtered sample_profile\n");
	      return FAILED;
	    }
	    if ((fp_stats = fopen(sim->profile_stats, "w+")) == NULL) {
	      fprintf(stderr, "ERROR: Cannot open stats sample_profile\n");
	      return FAILED;
	    }
      	break;
    case PROCESS_SAMPLING_REUSE:
	    if ((fp_filtered = fopen(sim->profile_fq, "r")) == NULL) {
	      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
	      return FAILED;
	    }
	    if ((fp_stats = fopen(sim->profile_stats, "r")) == NULL) {
	      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
	      return FAILED;
	    }
      	break;
    case PROCESS_MODEL:
    	if(set_model_qc(sim) == FAILED)
    	  return FAILED;

    	return SUCCEEDED;
      	break;
  }

  // these happen only if we aren't using model-basedc
  if (get_fastq_inf(sim, qc) == FAILED) {
    return FAILED;
  }
  print_fastq_stats(sim, fastq);

  return SUCCEEDED;
}

///////////////////////////////////////////////////////
// Function: parse_options - parse the user options  //
///////////////////////////////////////////////////////

void parse_options(int argc, char** argv, sim_t *sim)
{

  char *tp, *tmp_buf;
  long len;
  long num;
  long ratio;
  unsigned int seed = (unsigned int)time(NULL);

  // Variables for Option
  int opt, option_index;
  struct option long_options[] = {
    {"sample-fastq", 1, NULL, 0},
    {"data-type", 1, NULL, 0},
    {"depth", 1, NULL, 0},
    {"length-mean", 1, NULL, 0},
    {"length-sd", 1, NULL, 0},
    {"length-min", 1, NULL, 0},
    {"length-max", 1, NULL, 0},
    {"accuracy-mean", 1, NULL, 0},
    {"accuracy-sd", 1, NULL, 0},
    {"accuracy-min", 1, NULL, 0},
    {"accuracy-max", 1, NULL, 0},
    {"difference-ratio", 1, NULL, 0},
    {"model_qc", 1, NULL, 0},
    {"prefix", 1, NULL, 0},
    {"sample-profile-id", 1, NULL, 0},
    {"seed", 1, NULL, 0},
    {"help", 0, NULL, 0},
    {"version", 0, NULL, 0},
    {0, 0, 0, 0}
  };

  // Option parsing
  option_index = 0;
  while ((opt = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (opt) {
    case 0:
      sim->set_flg[option_index] = 1;

      switch (option_index) {
      case 0:
        if ((fastq.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(fastq.file, optarg);
        break;

      case 1:
        if (strncmp(optarg, "CLR", 3) == 0) {
          sim->data_type = DATA_TYPE_CLR;
        } else if (strncmp(optarg, "CCS", 3) == 0) {
          sim->data_type = DATA_TYPE_CCS;
        } else {
          fprintf(stderr, "ERROR (data-type: %s): Acceptable value is CLR or CCS.\n", optarg);
          exit(-1);
        }
        break;

      case 2:
        sim->depth = atof(optarg);
        if (sim->depth <= 0.0) {
          fprintf(stderr, "ERROR (depth: %s): Acceptable range is more than 0.\n", optarg);
          exit(-1);
        }
        break;

      case 3:
        sim->len_mean = atof(optarg);
        if ((sim->len_mean < 1) || (sim->len_mean > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-mean: %s): Acceptable range is 1-%ld.\n", optarg, (long)(FASTQ_LEN_MAX));
          exit(-1);
        }
        break;

      case 4:
        sim->len_sd = atof(optarg);
        if ((sim->len_sd < 0) || (sim->len_sd > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-sd: %s): Acceptable range is 0-%ld.\n",
            optarg, (long)(FASTQ_LEN_MAX));
          exit(-1);
        }
        break;

      case 5:
        if (strlen(optarg) >= 8) {
          fprintf(stderr, "ERROR (length-min: %s): Acceptable range is 1-%ld.\n", optarg, (long)(FASTQ_LEN_MAX));
          exit(-1);
        }
        sim->len_min = atoi(optarg);
        if ((sim->len_min < 1) || (sim->len_min > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-min: %s): Acceptable range is 1-%ld.\n", optarg, (long)(FASTQ_LEN_MAX));
          exit(-1);
        }
        break;

      case 6:
        if (strlen(optarg) >= 8) {
          fprintf(stderr, "ERROR (length-max: %s): Acceptable range is 1-%ld.\n", optarg, (long)(FASTQ_LEN_MAX));
          exit(-1);
        }
        sim->len_max = atoi(optarg);
        if ((sim->len_max < 1) || (sim->len_max > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-max: %s): Acceptable range is 1-%ld.\n", optarg, (long)(FASTQ_LEN_MAX));
          exit(-1);
        }
        break;

      case 7:
        sim->accuracy_mean = atof(optarg);
        if ((sim->accuracy_mean < 0.0) || (sim->accuracy_mean > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-mean: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 8:
        sim->accuracy_sd = atof(optarg);
        if ((sim->accuracy_sd < 0.0) || (sim->accuracy_sd > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-sd: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 9:
        sim->accuracy_min = atof(optarg);
        if ((sim->accuracy_min < 0.0) || (sim->accuracy_min > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-min: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 10:
        sim->accuracy_max = atof(optarg);
        if ((sim->accuracy_max < 0.0) || (sim->accuracy_max > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-max: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 11:
        if ((tmp_buf = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(tmp_buf, optarg);
        num = 0;
        tp = strtok(tmp_buf, ":");
        while (num < 3) {
          if (tp == NULL) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Format is sub:ins:del.\n", optarg);
            exit(-1);
          }
          if (strlen(tp) >= 5) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Acceptable range is 0-%d.\n", optarg, RATIO_MAX);
            exit(-1);
          }
          ratio = atoi(tp);
          if ((ratio < 0) || (ratio > RATIO_MAX)) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Acceptable range is 0-%d.\n", optarg, RATIO_MAX);
            exit(-1);
          }
          if (num == 0) {
            sim->sub_ratio = ratio;
          } else if (num == 1) {
            sim->ins_ratio = ratio;
          } else if (num == 2) {
            sim->del_ratio = ratio;
          }
          num ++;
          tp = strtok(NULL, ":");
        }
        free(tmp_buf);
        break;

      case 12:
        if ((sim->model_qc_file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim->model_qc_file, optarg);
        break;

      case 13:
        if ((sim->prefix = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim->prefix, optarg);
        break;

      case 14:
        if ((sim->profile_id = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim->profile_id, optarg);
        break;

      case 15:
        seed = (unsigned int)atoi(optarg);
        break;

      case 16: // help
        print_help();
        exit(0);

      case 17: // version
        printf("pbsim %s\n", versionString);
        exit(0);

      defalut:
        break;
      }

    default:
      break;
    }
  }

  srand((unsigned int)seed);
}

///////////////////////////////////////////////////////
// Function: get_ref_inf - Get reference information //
///////////////////////////////////////////////////////

int get_ref_inf(const sim_t *sim, ref_t *ref) {
  FILE *fp;
  char line[BUF_SIZE];
  long max_len = 0;
  int ret;

  fprintf(stderr, ":::: Reference stats ::::\n\n");
  fprintf(stderr, "file name : %s\n", ref->file);
  fprintf(stderr, "\n");

  if ((fp = fopen(ref->file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open reference file: %s\n", ref->file);
    return FAILED;
  }

  ref->num_seq = 0;
  ref->len = 0;

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    ret = trim(line);

    if (line[0] == '>') {
      if (ref->num_seq != 0) {
        if (ref->len < REF_SEQ_LEN_MIN) {
          fprintf(stderr, "ERROR: Reference is too short. Acceptable length >= %ld.\n", (long)(REF_SEQ_LEN_MIN));
          fclose(fp);
          return FAILED;
        }
        fprintf(stderr, "ref.%ld (len:%ld) : %s\n", ref->num_seq, ref->len, ref->id);
        fclose(fp_ref);
        if (ref->len > max_len) {
          max_len = ref->len;
        }
      }

      ref->num_seq ++;
      if (ref->num_seq > REF_SEQ_NUM_MAX) {
        fprintf(stderr, "ERROR: References are too many. Max number of reference is %ld.\n", (long)(REF_SEQ_NUM_MAX));
        fclose(fp);
        return FAILED;
      }

      strncpy(ref->id, &line[1], REF_ID_LEN_MAX);
      ref->id[std::min(REF_ID_LEN_MAX, (int)strlen(&line[1]))] = '\0';

      sprintf(sim->outfile_ref, "%s_%04ld.ref", sim->prefix, ref->num_seq);
      if ((fp_ref = fopen(sim->outfile_ref, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim->outfile_ref);
        fclose(fp);
        return FAILED;
      }

      ref->len = 0;

      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }

      fprintf(fp_ref, ">%s\n", ref->id);
    } else {
      ref->len += strlen(line);

      if (ref->len > REF_SEQ_LEN_MAX) {
        fprintf(stderr, "ERROR: Reference is too long. Acceptable length <= %ld.\n", (long)(REF_SEQ_LEN_MAX));
        fclose(fp);
        return FAILED;
      }

      fprintf(fp_ref, "%s\n", line);
    }
  }
  fclose(fp);

  if (ref->len < REF_SEQ_LEN_MIN) {
    fprintf(stderr, "ERROR: Reference is too short. Acceptable length >= %ld.\n", (long)(REF_SEQ_LEN_MIN));
    return FAILED;
  }
  fprintf(stderr, "ref.%ld (len:%ld) : %s\n", ref->num_seq, ref->len, ref->id);
  fclose(fp_ref);
  if (ref->len > max_len) {
    max_len = ref->len;
  }

  fprintf(stderr, "\n");

  if ((ref->seq = (char *)malloc(max_len + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////
// Function: get_ref_seq - Get reference sequence //
////////////////////////////////////////////////////

int get_ref_seq(const sim_t *sim, ref_t *ref) {
  FILE *fp;
  char line[BUF_SIZE];
  long offset = 0;
  long copy_size;
  int ret;

  sprintf(sim->outfile_ref, "%s_%04ld.ref", sim->prefix, ref->num);

  if ((fp = fopen(sim->outfile_ref, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open outfile_ref file: %s\n", sim->outfile_ref);
    return FAILED;
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    ret = trim(line);

    if (line[0] == '>') {
      memcpy(ref->id, &line[1], strlen(line)-1);
      ref->id[std::min(REF_ID_LEN_MAX, (int)strlen(&line[1]))] = '\0';
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }
    } else {
      copy_size = strlen(line);
      memcpy(ref->seq + offset, line, copy_size);
      offset += copy_size;
    }
  }
  fclose(fp);

  ref->seq[offset] = '\0';
  ref->len = strlen(ref->seq);

  return SUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: get_fastq_inf - Get FASTQ information //
/////////////////////////////////////////////////////

int get_fastq_inf(const sim_t *sim, const qc_t *qc) {
  FILE *fp;
  char *tp, *item;
  char line[BUF_SIZE];
  char qc_tmp[FASTQ_LEN_MAX];
  long len;
  double prob;
  double accuracy;
  double accuracy_total = 0;
  long value;
  double variance;
  long i;
  int line_num;

  for (i=0; i<=sim->len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  fastq.num = 0;
  fastq.len_min = LONG_MAX;
  fastq.len_max = 0;
  fastq.len_total = 0;
  fastq.num_filtered = 0;
  fastq.len_min_filtered = LONG_MAX;
  fastq.len_max_filtered = 0;
  fastq.len_total_filtered = 0;

  if (sim->process == PROCESS_SAMPLING_REUSE) {
    while (fgets(line, BUF_SIZE, fp_stats) != NULL) {
      trim(line);
      tp = strtok(line, "\t");
      item = tp;
      tp = strtok(NULL, "\t");

      if (strcmp(item, "num") == 0) {
        fastq.num_filtered = atol(tp);
      } else if (strcmp(item, "len_total") == 0) {
        fastq.len_total_filtered = atol(tp);
      } else if (strcmp(item, "len_min") == 0) {
        fastq.len_min_filtered = atol(tp);
      } else if (strcmp(item, "len_max") == 0) {
        fastq.len_max_filtered = atol(tp);
      } else if (strcmp(item, "len_mean") == 0) {
        fastq.len_mean_filtered = atof(tp);
      } else if (strcmp(item, "len_sd") == 0) {
        fastq.len_sd_filtered = atof(tp);
      } else if (strcmp(item, "accuracy_mean") == 0) {
        fastq.accuracy_mean_filtered = atof(tp);
      } else if (strcmp(item, "accuracy_sd") == 0) {
        fastq.accuracy_sd_filtered = atof(tp);
      }
    }
  } else {
    if ((fp = fopen(fastq.file, "r")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open fastq file: %s\n", fastq.file);
      return FAILED;
    }

    qc_tmp[0] = '\0';
    len = 0;
    line_num = 0;

    while (fgets(line, BUF_SIZE, fp) != NULL) {
      if (trim(line) == EXISTS_LINE_FEED) {
        line_num ++;

        if (line_num == 4) {
          len += strlen(line);

          if (len > FASTQ_LEN_MAX) {
            fprintf(stderr, "ERROR: fastq is too long. Max acceptable length is %ld.\n", (long)(FASTQ_LEN_MAX));
            fclose(fp);
            return FAILED;
          }

          fastq.num ++;
          fastq.len_total += len;

          if (fastq.num > FASTQ_NUM_MAX) {
            fprintf(stderr, "ERROR: fastq is too many. Max acceptable number is %ld.\n", (long)(FASTQ_NUM_MAX));
            fclose(fp);
            return FAILED;
          }

          if (len > fastq.len_max) {
            fastq.len_max = len;
          }
          if (len < fastq.len_min) {
            fastq.len_min = len;
          }

          if ((len >= sim->len_min) && (len <= sim->len_max)) {
            strcat(qc_tmp, line);
            prob = 0.0;
            for (i=0; i<len; i++) {
              prob += qc[(int)qc_tmp[i] - 33].prob;
            }
            accuracy = 1.0 - (prob / len);

            if ((accuracy >= sim->accuracy_min) && (accuracy <= sim->accuracy_max)) {
              accuracy_total += accuracy;
              fastq.num_filtered ++;
              fastq.len_total_filtered += len;

              freq_len[len] ++;
              value = (int)(accuracy * 100000 + 0.5);
              freq_accuracy[value] ++;

              fprintf(fp_filtered, "%s\n", qc_tmp);

              if (len > fastq.len_max_filtered) {
                fastq.len_max_filtered = len;
              }
              if (len < fastq.len_min_filtered) {
                fastq.len_min_filtered = len;
              }
            }
          }

          line_num = 0;
          qc_tmp[0] = '\0';
          len = 0;
        }
      } else {
        if (line_num == 3) {
          len += strlen(line);
          if (len > FASTQ_LEN_MAX) {
            fprintf(stderr, "ERROR: fastq is too long. Max acceptable length is %ld.\n", (long)(FASTQ_LEN_MAX));
            fclose(fp);
            return FAILED;
          }
          strcat(qc_tmp, line);
        }
      }
    }

    fclose(fp);

    if (fastq.num_filtered < 1) {
      fprintf(stderr, "ERROR: there is no sample-fastq in the valid range of length and accuracy.\n");
      return FAILED;
    }

    fastq.len_mean_filtered = (double)fastq.len_total_filtered / fastq.num_filtered;
    fastq.accuracy_mean_filtered = accuracy_total / fastq.num_filtered;

    variance = 0.0;
    for (i=0; i<=sim->len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((fastq.len_mean_filtered - i), 2) * freq_len[i];
      }
    }
    fastq.len_sd_filtered = sqrt(variance / fastq.num_filtered);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((fastq.accuracy_mean_filtered - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    fastq.accuracy_sd_filtered = sqrt(variance / fastq.num_filtered);

    if (sim->process == PROCESS_SAMPLING_STORE) {
      fprintf(fp_stats, "num\t%ld\n", fastq.num_filtered);
      fprintf(fp_stats, "len_total\t%lld\n", fastq.len_total_filtered);
      fprintf(fp_stats, "len_min\t%ld\n", fastq.len_min_filtered);
      fprintf(fp_stats, "len_max\t%ld\n", fastq.len_max_filtered);
      fprintf(fp_stats, "len_mean\t%f\n", fastq.len_mean_filtered);
      fprintf(fp_stats, "len_sd\t%f\n", fastq.len_sd_filtered);
      fprintf(fp_stats, "accuracy_mean\t%f\n", fastq.accuracy_mean_filtered);
      fprintf(fp_stats, "accuracy_sd\t%f\n", fastq.accuracy_sd_filtered);
    }
  }

  return SUCCEEDED;
}

//////////////////////////////////////////////////////////
// Function: init_sim_res - Initiate simulation results //
//////////////////////////////////////////////////////////

void init_sim_res(sim_t *sim) {
  sim->res_num = 0;
  sim->res_len_total = 0;
  sim->res_sub_num = 0;
  sim->res_ins_num = 0;
  sim->res_del_num = 0;
  sim->res_len_min = LONG_MAX;
  sim->res_len_max = 0;
}

/////////////////////////////////////////////////////////
// Function: set_sim_param - Set simulation parameters //
/////////////////////////////////////////////////////////

int set_sim_param(sim_t *sim) {
  FILE *fp;
  long sum;

  // data-type
  if (!(sim->set_flg[1])) {
    sim->data_type = DATA_TYPE_CLR;
  }

  // depth
  if (!(sim->set_flg[2])) {
    sim->depth = (sim->data_type == DATA_TYPE_CLR) ? 20.0 : 50.0;
  }

  // length-mean
  if (!(sim->set_flg[3])) {
    sim->len_mean = (sim->data_type == DATA_TYPE_CLR) ? 3000 : 450;
  }

  // length-sd
  if (!(sim->set_flg[4])) {
    sim->len_sd = (sim->data_type == DATA_TYPE_CLR) ? 2300 : 170;
  }

  // length-min
  if (!(sim->set_flg[5])) {
    sim->len_min = (sim->data_type == DATA_TYPE_CLR) ? 100 : 100;
  }

  // length-max
  if (!(sim->set_flg[6])) {
    sim->len_max = (sim->data_type == DATA_TYPE_CLR) ? 25000 : 2500;
  }

  // accuracy-mean
  if (sim->data_type == DATA_TYPE_CLR) {
    if (sim->set_flg[7]) {
      sim->accuracy_mean = int(sim->accuracy_mean * 100) * 0.01;
    } else {
      sim->accuracy_mean = 0.78;
    }
  } else {
    sim->accuracy_mean = 0.98;
  }

  // accuracy-sd
  if (sim->data_type == DATA_TYPE_CLR && sim->set_flg[8]) {
    sim->accuracy_sd = int(sim->accuracy_sd * 100) * 0.01;
  } else {
    sim->accuracy_sd = 0.02;
  }

  // accuracy-min
  if (sim->data_type == DATA_TYPE_CLR && sim->set_flg[9]) {
    sim->accuracy_min = int(sim->accuracy_min * 100) * 0.01;
  } else {
    sim->accuracy_min = 0.75;
  }

  // accuracy-max
  if (sim->data_type == DATA_TYPE_CLR && sim->set_flg[10]) {
    sim->accuracy_max = int(sim->accuracy_max * 100) * 0.01;
  } else {
    sim->accuracy_max = 1.0;
  }

  // difference-ratio
  if (!(sim->set_flg[11])) {
    if (sim->data_type == DATA_TYPE_CLR) {
      sim->sub_ratio = 10;
      sim->ins_ratio = 60;
      sim->del_ratio = 30;
    } else {
      sim->sub_ratio = 6;
      sim->ins_ratio = 21;
      sim->del_ratio = 73;
    }
  }

  sum = sim->sub_ratio + sim->ins_ratio + sim->del_ratio;
  sim->sub_rate = (double)sim->sub_ratio / sum;
  sim->ins_rate = (double)sim->ins_ratio / sum;
  sim->del_rate = (double)sim->del_ratio / sum;

  // prefix and outfile
  if (!(sim->set_flg[13])) {
    if ((sim->prefix = (char *)malloc(3)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      exit(-1);
    }
    strcpy(sim->prefix, "sd");
  }

  if ((sim->outfile_ref = (char *)malloc(strlen(sim->prefix) + 10)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((sim->outfile_fq = (char *)malloc(strlen(sim->prefix) + 12)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((sim->outfile_maf = (char *)malloc(strlen(sim->prefix) + 10)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  // profile
  if (sim->set_flg[14]) {
    if ((sim->profile_fq = (char *)malloc(strlen(sim->profile_id) + 22)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      return FAILED;
    }

    if ((sim->profile_stats = (char *)malloc(strlen(sim->profile_id) + 22)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      return FAILED;
    }

    sprintf(sim->profile_fq, "sample_profile_%s.fastq", sim->profile_id);
    sprintf(sim->profile_stats, "sample_profile_%s.stats", sim->profile_id);
  }

  // length and accuracy
  if (sim->len_min > sim->len_max) {
    fprintf(stderr, "ERROR: length min(%ld) is greater than max(%ld).\n", sim->len_min, sim->len_max);
    return FAILED;
  }
  if (sim->accuracy_min > sim->accuracy_max) {
    fprintf(stderr, "ERROR: accuracy min(%f) is greater than max(%f).\n", sim->accuracy_min, sim->accuracy_max);
    return FAILED;
  }

  // process
  if (sim->set_flg[12]) {
    if ((sim->set_flg[0]) || (sim->set_flg[14])) {
      fprintf(stderr, "ERROR: either --sample-fastq(and/or --sample-profile-id)(sampling-based) or --model_qc(model-based) should be set.\n");
      return FAILED;
    }
    sim->process = PROCESS_MODEL;
  } else {
    if (sim->set_flg[0]) {
      if (sim->set_flg[14]) {
        sim->process = PROCESS_SAMPLING_STORE;
      } else {
        sim->process = PROCESS_SAMPLING;
      }
    } else {
      if (sim->set_flg[14]) {
        sim->process = PROCESS_SAMPLING_REUSE;
      } else {
        fprintf(stderr, "ERROR: either --sample-fastq(and/or --sample-profile-id)(sampling-based) or --model_qc(model-based) should be set.\n");
        return FAILED;
      }
    }
  }

  // sample profile
  if (sim->process == PROCESS_SAMPLING_STORE) {
    if ((fp = fopen(sim->profile_fq, "r")) != NULL) {
      fprintf(stderr, "ERROR: %s exists.\n", sim->profile_fq);
      fclose(fp);
      return FAILED;
    }
    if ((fp = fopen(sim->profile_stats, "r")) != NULL) {
      fprintf(stderr, "ERROR: %s exists.\n", sim->profile_stats);
      fclose(fp);
      return FAILED;
    }
  }

  if (sim->process == PROCESS_SAMPLING_REUSE) {
    if ((fp = fopen(sim->profile_fq, "r")) == NULL) {
      fprintf(stderr, "ERROR: %s does not exist.\n", sim->profile_fq);
      return FAILED;
    }
    fclose(fp);
    if ((fp = fopen(sim->profile_stats, "r")) == NULL) {
      fprintf(stderr, "ERROR: %s does not exist.\n", sim->profile_stats);
      return FAILED;
    }
    fclose(fp);
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////////
// Function: simulate_by_sampling - Simulate by model //
////////////////////////////////////////////////////////

int simulate_by_sampling(sim_t *sim, ref_t *ref, mut_t *mut, fastq_t *fastq, qc_t *qc) {
  long len;
  long long len_total = 0;
  long sampling_num, sampling_interval, sampling_value, sampling_residue;
  long num;
  long i, j;
  long index;
  long value;
  double accuracy, accuracy_total = 0.0;
  double prob, variance;
  char id[128];
  int digit_num1[4], digit_num2[4], digit_num[4];

  for (i=0; i<=sim->len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  for (i=0; i<=93; i++) {
    mut->sub_thre[i] = int((qc[i].prob * sim->sub_rate) * 1000000 + 0.5);
    mut->ins_thre[i] = int((qc[i].prob * (sim->sub_rate + sim->ins_rate)) * 1000000 + 0.5);
  }
  mut->del_thre = int((1.0 - fastq->accuracy_mean_filtered) * sim->del_rate * 1000000 + 0.5);

  sampling_num = (long)(sim->len_quota / fastq->len_total_filtered);
  sampling_residue = sim->len_quota % fastq->len_total_filtered;
  if (sampling_residue == 0) {
    sampling_interval = 1;
  } else {
    sampling_interval = (long)((double)(fastq->len_total_filtered / sampling_residue) * 2 + 0.5);
    if (sampling_interval > (long)(fastq->num_filtered * 0.5)) {
      sampling_interval = (long)(fastq->num_filtered * 0.5);
    }
  }

  // Make simulation data
  while (len_total < sim->len_quota) {
    rewind(fp_filtered);

    sampling_value = rand() % fastq->num_filtered;
    while (fgets(mut->qc, fastq->len_max_filtered + 2, fp_filtered) != NULL) {
      if (len_total >= sim->len_quota) {
        break;
      }

      trim(mut->qc);

      if (sampling_value % sampling_interval == 0) {
        num = sampling_num + 1;
      } else {
        num = sampling_num;
      }
      sampling_value ++;

      for (i=0; i<num; i++) {
        if (len_total >= sim->len_quota) {
          break;
        }

        mut->tmp_len_max = sim->len_quota - len_total;
        if (mut->tmp_len_max < sim->len_min) {
          mut->tmp_len_max = sim->len_min;
        }

        if (mutate(sim, ref) == FAILED) {
          return FAILED;
        }

        sim->res_num ++;
        len = strlen(mut->new_seq);
        sim->res_len_total += len;
        len_total += len;
        freq_len[len] ++;

        if (len > sim->res_len_max) {
          sim->res_len_max = len;
        }
        if (len < sim->res_len_min) {
          sim->res_len_min = len;
        }

        prob = 0.0;
        for (j=0; j<len; j++) {
          prob += qc[(int)mut->new_qc[j] - 33].prob;
        }
        accuracy = 1.0 - (prob / len);
        accuracy_total += accuracy;
        value = (int)(accuracy * 100000 + 0.5);
        freq_accuracy[value] ++;

        sprintf(id, "S%ld_%ld", ref->num, sim->res_num);
        fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut->new_seq, id, mut->new_qc);

        digit_num1[0] = 3;
        digit_num2[0] = 1 + count_digit(sim->res_num);
        digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

        digit_num1[1] = count_digit((mut->seq_left - 1));
        digit_num2[1] = 1;
        digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

        digit_num1[2] = count_digit((mut->seq_right - mut->seq_left + 1));
        digit_num2[2] = count_digit(len);
        digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

        digit_num1[3] = count_digit(ref->len);
        digit_num2[3] = count_digit(len);
        digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

        fprintf(fp_maf, "a\ns ref");
        while (digit_num1[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num1[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld", mut->seq_left - 1);
        while (digit_num1[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld +", mut->seq_right - mut->seq_left + 1);
        while (digit_num1[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n", ref->len, mut->maf_ref_seq);
        fprintf(fp_maf, "s %s", id);
        while (digit_num2[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num2[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %d", 0);
        while (digit_num2[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %c", len, mut->seq_strand);
        while (digit_num2[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n\n", len, mut->maf_seq);
      }
    }

    sampling_num = 0;
  }

  sim->res_len_mean = (double)sim->res_len_total / sim->res_num;
  sim->res_accuracy_mean = accuracy_total / sim->res_num;

  if (sim->res_num == 1) {
    sim->res_len_sd = 0.0;
    sim->res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim->len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim->res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim->res_len_sd = sqrt(variance / sim->res_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim->res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim->res_accuracy_sd = sqrt(variance / sim->res_num);
  }

  return SUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: simulate_by_model - Simulate by Model //
/////////////////////////////////////////////////////

int simulate_by_model(sim_t *sim, ref_t *ref, qc_t *qc) {
  long len;
  long long len_total = 0;
  long num;
  long i, j, k;
  double prob, mean, variance, sd;
  double len_prob_total, accuracy_prob_total, qc_prob_total, value, sum;
  double accuracy_total = 0.0;
  int accuracy;
  long prob2len[100001], prob2accuracy[100001], prob2qc[101][1001];
  long len_rand_value, accuracy_rand_value, qc_rand_value[101];
  long start_wk, end_wk;
  long index;
  long accuracy_min, accuracy_max;
  char id[128];
  int digit_num1[4], digit_num2[4], digit_num[4];

  for (i=0; i<=sim->len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  for (i=0; i<=93; i++) {
    mut.sub_thre[i] = int((qc[i].prob * sim->sub_rate) * 1000000 + 0.5);
    mut.ins_thre[i] = int((qc[i].prob * (sim->sub_rate + sim->ins_rate)) * 1000000 + 0.5);
  }
  mut.del_thre = int((1.0 - sim->accuracy_mean) * sim->del_rate * 1000000 + 0.5);

  accuracy_min = (long)(sim->accuracy_min * 100);
  accuracy_max = (long)(sim->accuracy_max * 100);

  // length distribution
  variance = log(1 + pow((sim->len_sd / sim->len_mean) ,2));
  mean = log(sim->len_mean) - variance * 0.5;
  sd = sqrt(variance);

  if (sim->len_sd == 0.0) {
    prob2len[1] = int(sim->len_mean + 0.5);
    len_rand_value = 1;
  } else {
    start_wk = 1;
    len_prob_total = 0.0;
    for (i=sim->len_min; i<=sim->len_max; i++) {
      len_prob_total += exp(-1 * pow((log(i)-mean), 2) / 2 / variance) / sqrt(2*M_PI) / sd / i;
      end_wk = int(len_prob_total * 100000 + 0.5);
      if (end_wk > 100000) {
        end_wk = 100000;
      }

      for (j=start_wk; j<=end_wk; j++) {
        prob2len[j] = i;
      }

      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    len_rand_value = end_wk;
  }

  if (len_rand_value < 1) {
    fprintf(stderr, "ERROR: length parameters are not appropriate.\n");
    return FAILED;
  }

  // accuracy distribution
  if (sim->data_type == DATA_TYPE_CLR) {
    mean = sim->accuracy_mean * 100;
    sd = sim->accuracy_sd * 100;
    variance = pow(sd, 2);

    if (sd == 0.0) {
      prob2accuracy[1] = int(mean + 0.5);
      accuracy_rand_value = 1;
    } else {
      start_wk = 1;
      accuracy_prob_total = 0.0;
      for (i=accuracy_min; i<=accuracy_max; i++) {
        accuracy_prob_total += exp(-1 * pow(i - mean, 2) / 2 / variance) / sqrt(2 * M_PI) / sd;
        end_wk = int(accuracy_prob_total * 100000 + 0.5);
        if (end_wk > 100000) {
          end_wk = 100000;
        }

        for (j=start_wk; j<=end_wk; j++) {
          prob2accuracy[j] = i;
        }

        if (end_wk >= 100000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      accuracy_rand_value = end_wk;
    }
  } else {
    sum = 0;
    for (i=accuracy_min; i<=accuracy_max; i++) {
      sum += exp(0.5 * (i - 75));
    }

    start_wk = 1;
    accuracy_prob_total = 0.0;
    for (i=accuracy_min; i<=accuracy_max; i++) {
      accuracy_prob_total += exp(0.5 * (i - 75)) / sum;
      end_wk = int(accuracy_prob_total * 100000 + 0.5);
      if (end_wk > 100000) {
        end_wk = 100000;
      }

      for (j=start_wk; j<=end_wk; j++) {
        prob2accuracy[j] = i;
      }

      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    accuracy_rand_value = end_wk;
  }

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // quality code distributiin
  for (i=accuracy_min; i<=accuracy_max; i++) {
    start_wk = 1;
    qc_prob_total = 0.0;

    for (j=model_qc[i].min; j<=model_qc[i].max; j++) {
      qc_prob_total += model_qc[i].prob[j];
      end_wk = int(qc_prob_total * 1000 + 0.5);
      if (end_wk > 1000) {
        end_wk = 1000;
      }

      for (k=start_wk; k<=end_wk; k++) {
        prob2qc[i][k] = j;
      }

      if (end_wk >= 1000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    qc_rand_value[i] = end_wk;
  }

  // simulation
  while (len_total < sim->len_quota) {
    index = rand() % len_rand_value + 1;
    len = prob2len[index];
    if (len_total + len > sim->len_quota) {
      len = sim->len_quota - len_total;

      if (len < sim->len_min) {
        len = sim->len_min;
      }
    }

    mut.tmp_len_max = len;

    index = rand() % accuracy_rand_value + 1;
    accuracy = prob2accuracy[index];

    num = 0;
    while (num < len) {
      index = rand() % qc_rand_value[accuracy] + 1;
      index = prob2qc[accuracy][index];
      mut.qc[num ++] = qc[index].character;
      if (num >= len) {
        break;
      }
    }
    mut.qc[num] = '\0';

    if (mutate(sim, ref) == FAILED) {
      return FAILED;
    }

    len = strlen(mut.new_seq);
    sim->res_len_total += len;
    len_total += len;
    freq_len[len] ++;
    sim->res_num ++;

    if (len > sim->res_len_max) {
      sim->res_len_max = len;
    }
    if (len < sim->res_len_min) {
      sim->res_len_min = len;
    }

    prob = 0.0;
    for (i=0; i<len; i++) {
      prob += qc[(int)mut.new_qc[i] - 33].prob;
    }
    value = 1.0 - (prob / len);
    accuracy_total += value;
    accuracy = (int)(value * 100000 + 0.5);
    freq_accuracy[accuracy] ++;

    sprintf(id, "S%ld_%ld", ref->num, sim->res_num);
    fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.new_seq, id, mut.new_qc);

    digit_num1[0] = 3;
    digit_num2[0] = 1 + count_digit(sim->res_num);
    digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

    digit_num1[1] = count_digit((mut.seq_left - 1));
    digit_num2[1] = 1;
    digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

    digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
    digit_num2[2] = count_digit(len);
    digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

    digit_num1[3] = count_digit(ref->len);
    digit_num2[3] = count_digit(len);
    digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

    fprintf(fp_maf, "a\ns %s", ref->id);
    while (digit_num1[0] ++ < digit_num[0]) {
      fprintf(fp_maf, " ");
    }
    while (digit_num1[1] ++ < digit_num[1]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld", mut.seq_left - 1);
    while (digit_num1[2] ++ < digit_num[2]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
    while (digit_num1[3] ++ < digit_num[3]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %s\n", ref->len, mut.maf_ref_seq);
    fprintf(fp_maf, "s %s", id);
    while (digit_num2[0] ++ < digit_num[0]) {
      fprintf(fp_maf, " ");
    }
    while (digit_num2[1] ++ < digit_num[1]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %d", 0);
    while (digit_num2[2] ++ < digit_num[2]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
    while (digit_num2[3] ++ < digit_num[3]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
  }

  sim->res_len_mean = (double)sim->res_len_total / sim->res_num;
  sim->res_accuracy_mean = accuracy_total / sim->res_num;

  if (sim->res_num == 1) {
    sim->res_len_sd = 0.0;
    sim->res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim->len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim->res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim->res_len_sd = sqrt(variance / sim->res_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim->res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim->res_accuracy_sd = sqrt(variance / sim->res_num);
  }

  return SUCCEEDED;
}

/////////////////////////////////////////////////////////////////
// Function: set_mut - Set mutation parameters and varianeces  //
/////////////////////////////////////////////////////////////////

int set_mut(sim_t *sim) {
  mut.sub_nt_a = (char *)"TGC";
  mut.sub_nt_t = (char *)"AGC";
  mut.sub_nt_g = (char *)"ATC";
  mut.sub_nt_c = (char *)"ATG";
  mut.sub_nt_n = (char *)"ATGC";
  mut.ins_nt   = (char *)"ATGC";

  if ((mut.qc = (char *)malloc(sim->len_max + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.new_qc = (char *)malloc(sim->len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.tmp_qc = (char *)malloc(sim->len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.seq = (char *)malloc(sim->len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.new_seq = (char *)malloc(sim->len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.maf_seq = (char *)malloc(sim->len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.maf_ref_seq = (char *)malloc(sim->len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  return SUCCEEDED;
}

////////////////////////////////////
// Function: mutate - Mutate read //
////////////////////////////////////

int mutate(sim_t *sim, ref_t *ref) {
  char *line;
  char nt;
  long num;
  long i, j;
  long index;
  long rand_value;
  long qc_value;
  long len;
  long offset, seq_offset, maf_offset;

  len = strlen(mut.qc);
  if (mut.tmp_len_max < len) {
    len = mut.tmp_len_max;
  }

  // Place deletions
  offset = 0;
  for (i=0; i<len-1; i++) {
    mut.tmp_qc[offset ++] = mut.qc[i];
    if (rand() % 1000000 < mut.del_thre) {
      mut.tmp_qc[offset ++] = ' ';
      sim->res_del_num ++;
    }
  }
  mut.tmp_qc[offset ++] = mut.qc[len - 1];
  mut.tmp_qc[offset] = '\0';

  len = strlen(mut.tmp_qc);

  if (len >= ref->len) {
    offset = 0;
    len = ref->len;
  } else {
    offset = rand() % (ref->len - len + 1);
  }

  mut.seq_left = offset + 1;
  mut.seq_right = offset + len;

  if (1) { // sim->res_num % 2 == 0) {
    mut.seq_strand = '+';

    for (i=0; i<len; i++) {
      nt = toupper(ref->seq[offset + i]);
      mut.seq[i] = nt;
    }
  } else {
    mut.seq_strand = '-';

    for (i=0; i<len; i++) {
      nt = toupper(ref->seq[offset + i]);

      if (nt == 'A') {
        mut.seq[len-1-i] = 'T';
      } else if (nt == 'T') {
        mut.seq[len-1-i] = 'A';
      } else if (nt == 'G') {
        mut.seq[len-1-i] = 'C';
      } else if (nt == 'C') {
        mut.seq[len-1-i] = 'G';
      } else {
        mut.seq[len-1-i] = nt;
      }
    }
  }
  mut.seq[len] = '\0';

  // Place substitutions and insertions
  offset = 0;
  seq_offset = 0;
  maf_offset = 0;
  for (i=0; i<len; i++) {
    nt = mut.seq[seq_offset ++];

    if (mut.tmp_qc[i] == ' ') {
      mut.maf_seq[maf_offset] = '-';
      mut.maf_ref_seq[maf_offset] = nt;
      maf_offset ++;
      continue;
    }

    mut.new_qc[offset] = mut.tmp_qc[i];

    rand_value = rand() % 1000000;
    qc_value = (int)mut.tmp_qc[i] - 33;

    if (rand_value < mut.sub_thre[qc_value]) {
      sim->res_sub_num ++;
      index = rand() % 3;
      if (nt == 'A') {
        mut.new_seq[offset] = mut.sub_nt_a[index];
      } else if (nt == 'T') {
        mut.new_seq[offset] = mut.sub_nt_t[index];
      } else if (nt == 'G') {
        mut.new_seq[offset] = mut.sub_nt_g[index];
      } else if (nt == 'C') {
        mut.new_seq[offset] = mut.sub_nt_c[index];
      } else {
        index = rand() % 4;
        mut.new_seq[offset] = mut.sub_nt_n[index];
      }
      mut.maf_ref_seq[maf_offset] = nt;
    } else if (rand_value < mut.ins_thre[qc_value]) {
      sim->res_ins_num ++;
      index = rand() % 8;
      if (index >= 4) {
        mut.new_seq[offset] = nt;
      } else {
        mut.new_seq[offset] = mut.ins_nt[index];
      }
      seq_offset --;
      if (mut.seq_strand == '+') {
        mut.seq_right --;
      } else {
        mut.seq_left ++;
      }
      mut.maf_ref_seq[maf_offset] = '-';
    } else {
      mut.new_seq[offset] = nt;
      mut.maf_ref_seq[maf_offset] = nt;
    }
    mut.maf_seq[maf_offset] = mut.new_seq[offset];
    maf_offset ++;
    offset ++;
  }
  mut.new_qc[offset] = '\0';
  mut.new_seq[offset] = '\0';
  mut.maf_seq[maf_offset] = '\0';
  mut.maf_ref_seq[maf_offset] = '\0';

  if (mut.seq_strand == '-') {
    revcomp(mut.maf_seq);
    revcomp(mut.maf_ref_seq);
  }

  return SUCCEEDED;
}

///////////////////////////////////////////////////////
// Function: set_model_qc - Set quality code model   //
///////////////////////////////////////////////////////

int set_model_qc(const sim_t *sim) {
  FILE *fp;
  char line[BUF_SIZE];
  char *tp;
  long accuracy;
  int num;
  int i, j;

  if ((fp = fopen(sim->model_qc_file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open model_qc file: %s\n", sim->model_qc_file);
    return FAILED;
  }

  for (i=0; i<=100; i++) {
    for (j=0; j<=93; j++) {
      model_qc[i].prob[j] = 0.0;
    }
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    trim(line);

    tp = strtok(line, "\t");
    accuracy = atoi(tp);

    num = 0;
    tp = strtok(NULL, "\t");
    while (tp != NULL) {
      model_qc[accuracy].prob[num] = atof(tp);
      num ++;
      tp = strtok(NULL, "\t");
    }
  }
  fclose(fp);

  for (i=0; i<=100; i++) {
    model_qc[i].min = 0;
    model_qc[i].max = 93;

    for (j=0; j<=93; j++) {
      if (model_qc[i].prob[j] > 0.0) {
        model_qc[i].min = j;
        break;
      }
    }
    for (j=93; j>=0; j--) {
      if (model_qc[i].prob[j] > 0.0) {
        model_qc[i].max = j;
        break;
      }
    }
  }

  return SUCCEEDED;
}
