#ifndef PBWT_H
#define PBWT_H

#include "htslib/vcf.h"
#include "kvec.h"
#include "match.h"
#include "pbwt_col.h"
#include "phi.h"
#include "utils.h"
#include <htslib/synced_bcf_reader.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

typedef struct {
  uint32_t p_p;
  uint32_t p_l;
  uint32_t c_i;
  uint32_t c_r;
  uint32_t c_p;
  uint8_t c_s;
} q_state;

typedef struct {
  uint32_t n_haps;
  uint32_t n_sites;
  column_vec cols;
  phi phi;
  /* Mode-0-only: jump-distance threshold for choosing the backward LF-walk
   * over the Phi-based scan when computing a match length after a
   * mismatch (see pbwt_advance_site). Unused (never read) by the streamed
   * modes, which never have the full query panel in RAM. Set by
   * pbwt_query() from its thr90 argument; not populated by the streamed
   * query paths. */
  uint32_t thr90;
} pbwt;

void pbwt_build(char *filename, pbwt *pbwt, int threads);
void pbwt_update(uint8_t *col, uint32_t **pa, uint32_t **da,
                 uint32_t **spare_pa, uint32_t **spare_da, uint32_t n_h);

uv_res get_uv(const pbwt_col *col, uint32_t r);
static inline __attribute__((always_inline)) uint32_t get_r(const pbwt_col *col,
                                                            uint32_t i) {
  const c_arr *p = &col->p;
  uint32_t n = (uint32_t)p->n;
  if (n == 0)
    return 0;

  const uint8_t bits = p->bits;
  const uint32_t *data = p->data;

  if (i < c_arr_get(p, 0))
    return 0;
  if (i >= c_arr_get(p, n - 1))
    return n - 1;

  uint32_t pos = 0;
  uint32_t step = 1U << (31 - __builtin_clz(n));

  for (; step > 0; step >>= 1) {
    uint32_t next_pos = pos + step;
    if (next_pos < n) {
      uint64_t bit_off = (uint64_t)next_pos * bits;
      uint32_t idx = bit_off >> 5;
      uint32_t shift = bit_off & 31;

      uint32_t val = (data[idx] >> shift);
      if (shift + bits > 32) {
        val |= (data[idx + 1] << (32 - shift));
      }
      val &= ((1U << bits) - 1);

      if (val <= i) {
        pos = next_pos;
      }
    }
  }
  return pos;
}
static inline __attribute__((always_inline)) uint32_t
get_r_with_hint(const pbwt_col *col, uint32_t i, uint32_t hint) {
  const c_arr *p = &col->p;
  uint32_t n = (uint32_t)p->n;
  if (n == 0)
    return 0;

  uint32_t cur_val = c_arr_get(p, hint);
  if (cur_val <= i) {
    if (hint + 1 < n && c_arr_get(p, hint + 1) > i)
      return hint;
  } else {
    if (hint > 0 && c_arr_get(p, hint - 1) <= i)
      return hint - 1;
  }

  return get_r(col, i);
}
uint32_t fl(const pbwt_col *col, uint32_t i, uint32_t r);
uint32_t lf(const column_vec *cols, uint32_t col_idx, uint32_t i);

uint32_t get_l(const pbwt *pbwt, const c_arr *q, size_t q_base, uint32_t i,
              uint32_t c_i);
uint32_t get_l_u(const pbwt *pbwt, uint32_t p, uint32_t n, uint32_t c);
uint32_t get_l_d(const pbwt *pbwt, uint32_t p, uint32_t n, uint32_t c);

match_vec compute_smem(const pbwt *pbwt, const c_arr *q, size_t q_base,
                       uint32_t q_n);

int_vec get_haps_ms(const pbwt *pbwt, uint32_t p, uint32_t l, uint32_t c);
/* Empirically chosen (see thr90_sweep tuning experiment: a 9-panel x
 * 9-threshold sweep across panel_scaling/query_scaling/scale_stress,
 * cross-checked with a contention-free sequential spot-check) to minimize
 * worst-case regret vs. the per-panel optimum across widely varying panel
 * geometries. Still just a single fixed constant -- a strong correlation
 * with site count (r~=0.98 at fixed haplotype count) suggests a
 * size-adaptive default could do better, but that needs a dedicated
 * follow-up study before committing to a formula. */
#define MUPBWT_THR90_DEFAULT 60

/* thr90 < 0 means "use MUPBWT_THR90_DEFAULT"; only meaningful for mode 0
 * (the streamed query modes never use this threshold at all). */
void pbwt_query(const char *pbwt_filename, const char *filename, int threads,
                int thr90);
void pbwt_query_m(const char *pbwt_filename, const char *filename, int threads);
void pbwt_query_q(const char *pbwt_filename, const char *filename, int threads);
void pbwt_query_mq(const char *pbwt_filename, const char *filename,
                   int threads);
void free_pbwt(pbwt *p);
int pbwt_serialize(FILE *fp, const pbwt *p);
int pbwt_deserialize(FILE *fp, pbwt *p);
int pbwt_serialize_gz(const char *path, const pbwt *p);
int pbwt_deserialize_gz(const char *path, pbwt *p);
void print_pbwt(const pbwt *p);
void pbwt_opt(const char *filename, pbwt *pbwt);
void print_pa(const pbwt *p, uint32_t c);
void print_da(const pbwt *p, uint32_t c);
void print_back(const pbwt *pbwt, uint32_t c);
void print_for(const pbwt *pbwt, uint32_t c);
void print_lce(const pbwt *pbwt, uint32_t c);
void print_row(const pbwt *pbwt, uint32_t r);

static void pbwt_print_size(const pbwt *p) {
  if (!p)
    return;

  double cols_size_mb = 0;
  size_t n_r = 0;
  for (uint32_t i = 0; i < p->n_sites; ++i) {
    cols_size_mb += pbwt_col_size_mb(&p->cols.a[i]);
    n_r += p->cols.a[i].p.n;
  }

  double phi_size_mb_val = phi_size_mb(&p->phi);

  double overhead_mb = (double)sizeof(pbwt) / (1024.0 * 1024.0);

  double total_size_mb = cols_size_mb + phi_size_mb_val + overhead_mb;

  printf("PBWT Memory Usage Breakdown:\n");
  printf("----------------------------\n");
  printf("Number of Haplotypes: %u\n", p->n_haps);
  printf("Number of Sites:      %u\n", p->n_sites);
  printf("Number of Runs:       %u\n", n_r);
  printf("Columns (BWT/RLE):    %.4f MB\n", cols_size_mb);
  printf("Phi Structures:       %.4f MB\n", phi_size_mb_val);
  printf("Struct Overhead:      %.4f MB\n", overhead_mb);
  printf("----------------------------\n");
  printf("Total PBWT Size:      %.4f MB\n", total_size_mb);
}

#endif
