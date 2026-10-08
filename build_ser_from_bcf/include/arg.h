#ifndef ARG_H
#define ARG_H

#include <argp.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct arguments {
  char *input_file;
  char *output_file;
  char *query_file;
  int threads;
  int mode;
  int thr90;
};

static struct argp index_argp;
static struct argp query_argp;
static struct argp stats_argp;

static struct argp_option index_options[] = {
    {"input", 'i', "FILE", 0, "Input BCF file", 0},
    {"output", 'o', "FILE", 0, "Output index file", 0},
    {"threads", 't', "N", 0, "Number of threads", 0},
    {0}};

static error_t index_parser(int key, char *arg, struct argp_state *state) {
  struct arguments *args = state->input;
  switch (key) {
  case 'i':
    args->input_file = arg;
    break;
  case 'o':
    args->output_file = arg;
    break;
  case 't':
    args->threads = atoi(arg);
    if (args->threads <= 0) {
      args->threads = 1;
      fprintf(stderr, "Negative number of threads, threads setted to 1\n");
    } else if (args->threads > omp_get_max_threads()) {
      args->threads = omp_get_max_threads();
      fprintf(stderr, "Max number of threads = %d, threads setted to %d\n",
              omp_get_max_threads(), args->threads);
    }
    break;
  case ARGP_KEY_END:
    if (!args->input_file || !args->output_file) {
      fprintf(stderr, "Missing required options\n\n");
      argp_help(&index_argp, stdout, ARGP_HELP_STD_HELP, state->argv[0]);
      exit(1);
    }
    break;

  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp index_argp = {.options = index_options,
                                 .parser = index_parser,
                                 .args_doc = NULL,
                                 .doc = "index command",
                                 .children = NULL,
                                 .help_filter = NULL,
                                 .argp_domain = NULL};

static struct argp_option query_options[] = {
    {"input", 'i', "FILE", 0, "Input index file"},
    {"query", 'q', "FILE", 0, "Query BCF file"},
    {"threads", 't', "INT", 0, "Number of OpenMP threads (default 1)"},
    {"mode", 'm', "INT", 0,
     "Streaming mode: 0 = None (All in RAM), 1 = PBWT online (m), 2 = Query "
     "online (q), 3 = Both online (mq). Default: 0"},
    {"thr", 'T', "INT", 0,
     "Mode-0-only jump-distance threshold for choosing the backward LF-walk "
     "over the Phi-based scan after a mismatch. Ignored by modes 1-3. "
     "Default: 60"},
    {0}};

static error_t query_parser(int key, char *arg, struct argp_state *state) {
  struct arguments *args = state->input;
  switch (key) {
  case 'i':
    args->input_file = arg;
    break;
  case 'q':
    args->query_file = arg;
    break;
  case 't':
    args->threads = atoi(arg);
    if (args->threads <= 0) {
      args->threads = 1;
      fprintf(stderr, "Negative number of threads, threads set to 1\n");
    } else if (args->threads > omp_get_max_threads()) {
      args->threads = omp_get_max_threads();
      fprintf(stderr, "Max number of threads = %d, threads set to %d\n",
              omp_get_max_threads(), args->threads);
    }
    break;
  case 'm':
    args->mode = atoi(arg);
    if (args->mode < 0 || args->mode > 3) {
      fprintf(stderr, "Invalid mode %d. Must be 0, 1, 2, or 3.\n", args->mode);
      exit(EINVAL);
    }
    break;
  case 'T':
    args->thr90 = atoi(arg);
    if (args->thr90 < 0) {
      fprintf(stderr, "Invalid --thr %d: must be >= 0.\n", args->thr90);
      exit(EINVAL);
    }
    break;
  case ARGP_KEY_END:
    if (!args->input_file || !args->query_file) {
      fprintf(stderr, "Missing required options\n\n");
      argp_help(&query_argp, stdout, ARGP_HELP_STD_HELP, state->argv[0]);
      exit(1);
    }

    break;

  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp query_argp = {.options = query_options,
                                 .parser = query_parser,
                                 .args_doc = NULL,
                                 .doc = "query command",
                                 .children = NULL,
                                 .help_filter = NULL,
                                 .argp_domain = NULL};

static struct argp_option stats_options[] = {
    {"input", 'i', "FILE", 0, "Input index file", 0},
    {"threads", 't', "N", 0, "Number of threads", 0},
    {0}};

static error_t stats_parser(int key, char *arg, struct argp_state *state) {
  struct arguments *args = state->input;
  switch (key) {
  case 'i':
    args->input_file = arg;
    break;
  case 't':
    args->threads = atoi(arg);
    if (args->threads <= 0) {
      args->threads = 1;
      fprintf(stderr, "Negative number of threads, threads setted to 1\n");
    } else if (args->threads > omp_get_max_threads()) {
      args->threads = omp_get_max_threads();
      fprintf(stderr, "Max number of threads = %d, threads setted to %d\n",
              omp_get_max_threads(), args->threads);
    }
    break;

  case ARGP_KEY_END:
    if (!args->input_file) {
      fprintf(stderr, "Missing required options\n\n");
      argp_help(&stats_argp, stdout, ARGP_HELP_STD_HELP, state->argv[0]);
      exit(1);
    }
    break;

  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp stats_argp = {.options = stats_options,
                                 .parser = stats_parser,
                                 .args_doc = NULL,
                                 .doc = "stats command",
                                 .children = NULL,
                                 .help_filter = NULL,
                                 .argp_domain = NULL};
#endif
