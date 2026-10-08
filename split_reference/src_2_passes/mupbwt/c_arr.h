#ifndef C_ARR_H
#define C_ARR_H

#include "kvec.h"
#include "utils.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

typedef struct {
  uint32_t *data;
  size_t n;
  uint8_t bits;
} c_arr;

c_arr *c_arr_create(size_t n, uint32_t max_val);
c_arr *c_arr_from_array(const uint32_t *a, size_t n);
c_arr *c_arr_from_kvec(const int_vec kv);
c_arr *c_arr_from_array_max(const uint32_t *a, size_t n, uint32_t max_val);
c_arr *c_arr_from_kvec_max(const int_vec kv, uint32_t max_val);
void c_arr_free(c_arr *a);
void c_arr_set(c_arr *a, size_t i, uint32_t val);
// uint32_t c_arr_get(const c_arr *a, size_t i);
static inline __attribute__((always_inline)) uint32_t
c_arr_get(const c_arr *arr, size_t i) {
  if (i >= arr->n)
    return 0;
  size_t bit_pos = i * arr->bits;
  size_t idx = bit_pos / 32;
  size_t offset = bit_pos % 32;
  uint32_t val = (arr->data[idx] >> offset) & ((1u << arr->bits) - 1);
  if (offset + arr->bits > 32) {
    int bits_left = (offset + arr->bits) - 32;
    val |= (arr->data[idx + 1] & ((1u << bits_left) - 1))
           << (arr->bits - bits_left);
  }
  return val;
}

bool c_arr_equal(const c_arr *a, const c_arr *b);
void c_arr_print(const c_arr *a);
uint32_t c_arr_rank(const c_arr *a, uint32_t v);

static void print_c_arr_as_quad(const c_arr *a) {
  for (size_t i = 0; i + 3 < a->n; i += 4)
    printf("[%u %u %u %u] ", c_arr_get(a, i), c_arr_get(a, i + 1),
           c_arr_get(a, i + 2), c_arr_get(a, i + 3));
  printf("\n");
}

static void print_c_arr_as_penta(const c_arr *a) {
  for (size_t i = 0; i + 4 < a->n; i += 5)
    printf("[%u %u %u %u %u] ", c_arr_get(a, i), c_arr_get(a, i + 1),
           c_arr_get(a, i + 2), c_arr_get(a, i + 3), c_arr_get(a, i + 4));
  printf("\n");
}

int c_arr_serialize(const c_arr *arr, FILE *fp);
c_arr *c_arr_deserialize(FILE *fp);
int c_arr_serialize_gz(gzFile fp, const c_arr *arr);
int c_arr_deserialize_inplace_gz(gzFile fp, c_arr *arr);

static size_t c_arr_buffer_bytes(const c_arr *a) {
  if (a->n == 0 || a->bits == 0)
    return 0;
  uint64_t total_bits = (uint64_t)a->n * a->bits;
  uint64_t num_words = (total_bits + 31) / 32;
  return (size_t)num_words * sizeof(uint32_t);
}

static size_t c_arr_total_bytes(const c_arr *a) {
  if (!a)
    return 0;

  size_t header_size = sizeof(c_arr);

  uint64_t total_bits = (uint64_t)a->n * a->bits;
  uint64_t num_words = (total_bits + 31) / 32;
  size_t buffer_size = (size_t)num_words * sizeof(uint32_t);

  return header_size + buffer_size;
}

static double c_arr_size_mb(const c_arr *a) {
  if (!a)
    return 0.0;

  size_t struct_size = sizeof(c_arr);

  uint64_t total_bits = (uint64_t)a->n * a->bits;
  uint64_t num_words = (total_bits + 31) / 32;
  size_t buffer_size = (size_t)num_words * sizeof(uint32_t);

  size_t total_bytes = struct_size + buffer_size;

  return (double)total_bytes / (1024.0 * 1024.0);
}
#endif
