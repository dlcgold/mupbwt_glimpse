#include "c_arr.h"
#include "klib/kvec.h"
#include "pbwt.h"
#include "pbwt_col.h"
#include "phi.h"
#include "utils.h"
#include <arg.h>
#include <argp.h>
#include <stdio.h>
#include <sys/resource.h>
#include <time.h>

int main(int argc, char **argv) {

  if (argc < 2) {
    fprintf(stderr, "Usage: %s <index|query|stats> [options]\n", argv[0]);
    fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
    return 1;
  }

  struct arguments args = {0};
  char *cmd = argv[1];

  if (strcmp(cmd, "index") == 0) {
    argp_parse(&index_argp, argc - 1, argv + 1, 0, 0, &args);
    if (args.threads <= 0) {
      args.threads = 1;
      fprintf(stderr, "Threads not specified: using 1 thread\n");
    }
    // printf("INDEX: %s -> %s\n", args.input_file, args.output_file);
    fprintf(stderr, "mupbwt index with %d threads\n", args.threads);
    pbwt pbwt;
    struct timespec STARTM, ENDM;
    clock_gettime(CLOCK_MONOTONIC, &STARTM);
    pbwt_build(args.input_file, &pbwt, args.threads);
    clock_gettime(CLOCK_MONOTONIC, &ENDM);
    double wall_time =
        ENDM.tv_sec - STARTM.tv_sec + (ENDM.tv_nsec - STARTM.tv_nsec) / 1e9;

    fprintf(stderr, "mupbwt construction in %fs\n", wall_time);

    // print_pbwt(&pbwt);
    clock_t START = clock();
    FILE *out = fopen(args.output_file, "wb");
    if (!out) {
      perror("fopen");
      return 1;
    }
    if (pbwt_serialize(out, &pbwt) != 0) {
      fprintf(stderr, "Serialization error\n");
      fclose(out);
      return 1;
    }
    /* if (pbwt_serialize_gz(args.output_file, &pbwt) != 0) { */
    /*   fprintf(stderr, "Gzip Serialization error\n"); */
    /*   return 1; */
    /* } */
    /* fclose(out); */
    fprintf(stderr, "mupbwt disk dump in %fs\n",
            (float)(clock() - START) / CLOCKS_PER_SEC);
    free_pbwt(&pbwt);
  } else if (strcmp(cmd, "query") == 0) {
    args.mode = 0;
    args.thr90 = -1;

    argp_parse(&query_argp, argc - 1, argv + 1, ARGP_IN_ORDER, 0, &args);

    if (args.threads <= 0) {
      args.threads = 1;
      fprintf(stderr, "Threads not specified: using 1 thread\n");
    }

    fprintf(stderr, "mupbwt query with %d threads\n", args.threads);

    clock_t START = clock();
    struct timespec STARTM, ENDM;
    clock_gettime(CLOCK_MONOTONIC, &STARTM);

    switch (args.mode) {
    case 0:
      fprintf(stderr, "Run in mode 0: MuPBWT and Query in RAM\n");
      pbwt_query(args.input_file, args.query_file, args.threads, args.thr90);
      break;

    case 1:
      fprintf(stderr, "Run in mode 1: MuPBWT online, Query in RAM\n");
      pbwt_query_m(args.input_file, args.query_file, args.threads);
      break;

    case 2:
      fprintf(stderr, "Run in mode 2: MuPBWT in RAM, Query online\n");
      pbwt_query_q(args.input_file, args.query_file, args.threads);
      break;

    case 3:
      fprintf(stderr, "Run in mode 3: MuPBWT and Query online\n");
      pbwt_query_mq(args.input_file, args.query_file, args.threads);
      break;

    default:
      fprintf(stderr, "Unknown execution mode!\n");
      exit(-1);
    }

    clock_gettime(CLOCK_MONOTONIC, &ENDM);
    double wall_time =
        ENDM.tv_sec - STARTM.tv_sec + (ENDM.tv_nsec - STARTM.tv_nsec) / 1e9;

    fprintf(stderr, "Total time elapsed: %fs\n", wall_time);
  } else if (strcmp(cmd, "stats") == 0) {
    argp_parse(&stats_argp, argc - 1, argv + 1, 0, 0, &args);
    if (args.threads <= 0) {
      args.threads = 1;
      fprintf(stderr, "Threads not specified: using 1 thread\n");
    }

    fprintf(stderr, "mupbwt stats with %d threads\n", args.threads);
    clock_t START = clock();
    pbwt pbwt;
    /* FILE *in = fopen(args.input_file, "rb"); */
    /* pbwt_deserialize(in, &pbwt); */
    /* fclose(in); */
    if (pbwt_deserialize_gz(args.input_file, &pbwt) != 0) {
      fprintf(stderr, "load error %s\n", args.input_file);
      return 1;
    }

    fprintf(stderr, "mupbwt load from disk in %fs\n",
            (float)(clock() - START) / CLOCKS_PER_SEC);

    pbwt_print_size(&pbwt);
    free_pbwt(&pbwt);
  } else {
    fprintf(stderr, "Unknown command '%s'\n", cmd);
    fprintf(stderr, "Try '%s --help'\n", argv[0]);
    return 1;
  }

  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);

  fprintf(stderr, "Peak RAM usage: %.2f MB\n",
          (double)r_usage.ru_maxrss / 1024.0);
  return 0;
}
