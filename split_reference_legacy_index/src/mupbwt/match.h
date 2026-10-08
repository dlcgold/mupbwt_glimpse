#ifndef MATCH_H
#define MATCH_H

#include "utils.h"

typedef struct {
  uint32_t e;
  uint32_t l;
  uint32_t q_n;
  int_vec h;
} match;

typedef kvec_t(match) match_vec;

#define KV_PUSH_MATCH_VEC(vec, m)                                              \
  do {                                                                         \
    kv_push(match, (vec), (m));                                                \
  } while (0)

static inline void match_vec_free(match_vec *v) {
  for (size_t i = 0; i < v->n; i++) {
    kv_destroy(v->a[i].h);
  }
  kv_destroy(*v);
}

static inline void match_vec_print(const match_vec *v) {
  for (size_t i = 0; i < v->n; i++) {
    const match *m = &v->a[i];
    for (size_t j = 0; j < m->h.n; j++) {
      printf("MATCH\t%d\t%d\t%u\t%u\t%u\n", m->q_n, m->h.a[j],
             m->e - (m->l - 1), m->e, m->l);
    }
  }
}
#endif
