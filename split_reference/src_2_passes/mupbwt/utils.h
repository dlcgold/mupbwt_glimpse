#ifndef UTILS_H
#define UTILS_H

#include "kvec.h"
#include "stdlib.h"
#include "string.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#ifdef DEBUG
#define DBG(fmt, ...)                                                          \
  do {                                                                         \
    printf(fmt "\n", ##__VA_ARGS__);                                           \
  } while (0)
#else
#define DBG(fmt, ...)
#endif

#ifdef DEBUG
#define DBGL(fmt, ...)                                                         \
  do {                                                                         \
    printf(fmt " ", ##__VA_ARGS__);                                            \
  } while (0)
#else
#define DBGL(fmt, ...)
#endif

#define PRINT_ARRAY(fmt, arr, len)                                             \
  do {                                                                         \
    printf("[");                                                               \
    for (size_t _i = 0; _i < (len); _i++) {                                    \
      printf(fmt, (arr)[_i]);                                                  \
      if (_i + 1 < (len))                                                      \
        printf(", ");                                                          \
    }                                                                          \
    printf("]\n");                                                             \
  } while (0)

#ifdef DEBUG
#define DBG_PRINT_ARRAY(fmt, arr, len)                                         \
  do {                                                                         \
    DBGL("[");                                                                 \
    for (size_t _i = 0; _i < (len); _i++) {                                    \
      DBGL(fmt, (arr)[_i]);                                                    \
    }                                                                          \
    DBG("]\n");                                                                \
  } while (0)
#else
#define DBG_PRINT_ARRAY(fmt, arr, len)
#endif

#ifdef DEBUG
#define DBG_PRINT_KARRAY(fmt, arr)                                             \
  do {                                                                         \
    DBGL("[");                                                                 \
    for (size_t _i = 0; _i < (arr.n); _i++) {                                  \
      DBGL(fmt, (arr.a)[_i]);                                                  \
    }                                                                          \
    DBG("]");                                                                  \
  } while (0)
#else
#define DBG_PRINT_KARRAY(fmt, arr)
#endif

#define U32_MIN(a, b) ((uint32_t)((a) < (b) ? (a) : (b)))

typedef kvec_t(uint32_t) int_vec;
typedef kvec_t(uint8_t) int8_vec;
typedef kvec_t(int_vec) int_mat;

#define KV_PUSH_INT_VEC(mat, vec)                                              \
  do {                                                                         \
    kv_push(int_vec, mat, vec);                                                \
    kv_init(vec);                                                              \
  } while (0)

static void print_kvec_uint32_name(const int_vec *v, const char *name) {
  printf("%s: ", name);
  for (size_t i = 0; i < v->n; i++) {
    printf("%u ", v->a[i]);
  }
  printf("\n");
}

static void print_kvec_uint32(const int_vec *v) {
  for (size_t i = 0; i < v->n; i++) {
    printf("%u ", v->a[i]);
  }
  printf("\n");
}

static int cmp_uint32_t(const void *a, const void *b) {
  uint32_t ua = *(const uint32_t *)a;
  uint32_t ub = *(const uint32_t *)b;
  return (ua > ub) - (ua < ub);
}

static void sort_int_vec(int_vec *v) {
  if (v->n > 1) {
    qsort(v->a, v->n, sizeof(uint32_t), cmp_uint32_t);
  }
}

static inline __attribute__((always_inline)) uint8_t get_ns(bool z,
                                                            uint32_t r) {
  return (r & 1) ^ (!z);
}

typedef struct {
  uint32_t i;
  uint32_t r;
} map_res;

typedef struct {
  uint32_t r;
  uint32_t l;
} lce_res;

typedef struct {
  uint32_t u;
  uint32_t v;
} uv_res;

#endif
