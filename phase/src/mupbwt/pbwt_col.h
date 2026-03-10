
#ifndef PBWT_COL_H
#define PBWT_COL_H

#include "c_arr.h"
#include "kvec.h"
#include "zlib.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>

#define KV_PUSH_COL_VEC(vec, col)                                              \
  do {                                                                         \
    kv_push(pbwt_col, (vec), (col));                                           \
  } while (0)

typedef struct {
  bool zero;
  uint32_t n_zeros;

  c_arr p;
  c_arr uv;

  c_arr t;

  c_arr b_pa;
  c_arr b_da;
  c_arr e_pa;

} pbwt_col;

typedef kvec_t(pbwt_col) column_vec;

pbwt_col build_col(const uint8_t *c, const uint32_t *pa, const uint32_t *da,
                   uint32_t n_h);
void print_pbwt_col(const pbwt_col *col, uint32_t n_h);

void free_pbwt_col(pbwt_col *col);

int pbwt_col_serialize(FILE *fp, const pbwt_col *col);
int pbwt_col_deserialize(FILE *fp, pbwt_col *col);

int pbwt_col_serialize_gz(gzFile fp, const pbwt_col *col);
int pbwt_col_deserialize_gz(gzFile fp, pbwt_col *col);

static double pbwt_col_size_mb(const pbwt_col *col) {
  if (!col)
    return 0.0;

  size_t total_bytes = sizeof(pbwt_col);

  total_bytes += c_arr_buffer_bytes(&col->p);
  total_bytes += c_arr_buffer_bytes(&col->uv);
  total_bytes += c_arr_buffer_bytes(&col->t);
  total_bytes += c_arr_buffer_bytes(&col->b_pa);
  total_bytes += c_arr_buffer_bytes(&col->b_da);
  total_bytes += c_arr_buffer_bytes(&col->e_pa);

  return (double)total_bytes / (1024.0 * 1024.0);
}
#endif
