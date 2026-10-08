#ifndef PBWT_H
#define PBWT_H

#include "../containers/variant_map.h"
#include "../io/ref_genotype_reader.h"
#include "../utils/otools.h"
#include "../utils/verbose.h"
#include "htslib/vcf.h"
#include "kvec.h"
#include "match.h"
#include "pbwt_col.h"
#include "phi.h"
#include "utils.h"
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/regidx.h>
#include <htslib/synced_bcf_reader.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <zlib.h>

// #ifdef __cplusplus
// extern "C" {
// #endif

typedef struct {
  uint32_t p_p;
  uint32_t p_l;
  uint32_t c_i;
  uint32_t c_r;
  uint32_t c_p;
  uint8_t c_s;
  uint64_t *q_mask;
} smem_group;

typedef kvec_t(smem_group) group_vec;
typedef kvec_t(uint64_t) mask_arena;

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
  phi phid;
} pbwt;

void pbwt_build(char *filename, pbwt *pbwt, int threads);
void pbwt_build_af(pbwt *pbwt, std::string fref, ref_haplotype_set &H,
                   variant_map &V, std::string &region, bool common);
void pbwt_update(uint8_t *col, uint32_t **pa, uint32_t **da, uint32_t n_h);

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

uint32_t lf_r(const pbwt_col *col, uint32_t i, uint32_t r);
uint32_t get_l(const pbwt *pbwt, const uint8_t *q, uint32_t i, uint32_t c_i);
uint32_t get_l_pre(const pbwt *pbwt, const uint8_t *q, uint32_t i, uint32_t c_i,
                   uint32_t l);
uint32_t get_l_u(const pbwt *pbwt, uint32_t p, uint32_t n, uint32_t c);
uint32_t get_l_d(const pbwt *pbwt, uint32_t p, uint32_t n, uint32_t c);
uint32_t pbwt_lce(const pbwt *pbwt, uint32_t i1, uint32_t r1, uint32_t i2,
                  uint32_t r2, uint32_t c_start);
lce_res pbwt_lce_best(const pbwt *pbwt, uint32_t i1, uint32_t r1, uint32_t i2,
                      uint32_t r2, uint32_t i3, uint32_t r3, uint32_t c_start);

match_vec compute_smem(const pbwt *pbwt, const uint8_t *q, uint32_t q_n);

int_vec get_haps_ms(const pbwt *pbwt, uint32_t p, uint32_t l, uint32_t c);
void pbwt_query(const char *pbwt_filename, const char *filename, int threads);
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
  for (uint32_t i = 0; i < p->n_sites; ++i) {
    cols_size_mb += pbwt_col_size_mb(&p->cols.a[i]);
  }

  double phi_size_mb_val = phi_size_mb(&p->phid);

  double overhead_mb = (double)sizeof(pbwt) / (1024.0 * 1024.0);

  double total_size_mb = cols_size_mb + phi_size_mb_val + overhead_mb;

  printf("PBWT Memory Usage Breakdown:\n");
  printf("----------------------------\n");
  printf("Number of Haplotypes: %u\n", p->n_haps);
  printf("Number of Sites:      %u\n", p->n_sites);
  printf("Columns (BWT/RLE):    %.4f MB\n", cols_size_mb);
  printf("Phi Structures:       %.4f MB\n", phi_size_mb_val);
  printf("Struct Overhead:      %.4f MB\n", overhead_mb);
  printf("----------------------------\n");
  printf("Total PBWT Size:      %.4f MB\n", total_size_mb);
}

// #ifdef __cplusplus
// }
// #endif
#endif
