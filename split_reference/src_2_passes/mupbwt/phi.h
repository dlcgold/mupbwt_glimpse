#ifndef PHI_H
#define PHI_H

#include "c_arr.h"
#include "kvec.h"
#include "utils.h"

typedef struct {
  c_arr *phi_pos;
  c_arr *phi_inv_pos;
  c_arr *phi_supp;
  c_arr *phi_inv_supp;
  c_arr *phi_l_supp;
  uint32_t n_h;
  uint32_t n_w;
} phi;

void build_phi(phi *phi, int_vec *supp_b, int_vec *supp_e, int_vec *supp_pa_b,
               int_vec *supp_pa_e, int_vec *supp_da_b, uint32_t n_h,
               uint32_t n_w);

void phi_print(const phi *p);
void phi_free(phi *p);

int phi_serialize(const phi *p, FILE *fp);
int phi_deserialize(phi *p, FILE *fp);
int phi_serialize_gz(const phi *p, gzFile fp);
int phi_deserialize_gz(phi *p, gzFile fp);

// uint32_t phi_f(const phi *p, uint32_t pa_v, uint32_t c);
// uint32_t phi_inv_f(const phi *p, uint32_t pa_v, uint32_t c);
// uint32_t phi_l(const phi *p, uint32_t pa_v, uint32_t c);
static inline __attribute__((always_inline)) uint32_t phi_f(const phi *p,
                                                            uint32_t pa_v,
                                                            uint32_t c) {
  uint32_t t_c = c_arr_rank(&p->phi_pos[pa_v], c);
  t_c -= (t_c == p->phi_supp[pa_v].n) & (t_c > 0);
  return c_arr_get(&p->phi_supp[pa_v], t_c);
}

static inline __attribute__((always_inline)) uint32_t phi_inv_f(const phi *p,
                                                                uint32_t pa_v,
                                                                uint32_t c) {
  uint32_t t_c = c_arr_rank(&p->phi_inv_pos[pa_v], c);
  DBG("rank_%d(%d) = %d", c, pa_v, t_c);
  t_c -= (t_c == p->phi_inv_supp[pa_v].n) & (t_c > 0);
  DBG("update rank_%d(%d) = %d", c, pa_v, t_c);
  DBG("ret = %d", c_arr_get(&p->phi_inv_supp[pa_v], t_c));
  return c_arr_get(&p->phi_inv_supp[pa_v], t_c);
}

static inline __attribute__((always_inline)) uint32_t phi_l(const phi *p,
                                                            uint32_t pa_v,
                                                            uint32_t c) {
  if (c == 0 || phi_f(p, pa_v, c) == p->n_h)
    return 0;

  uint32_t t_c = c_arr_rank(&p->phi_pos[pa_v], c);
  t_c -= (t_c == p->phi_supp[pa_v].n) & (t_c > 0);

  uint32_t e_c = p->n_w - 1;
  e_c =
      e_c -
      ((t_c < p->phi_supp[pa_v].n - 1) ? e_c - c_arr_get(&p->phi_pos[pa_v], t_c)
                                       : 0) +
      ((t_c < p->phi_supp[pa_v].n - 1) ? 0 : e_c - e_c);

  return c_arr_get(&p->phi_l_supp[pa_v], t_c) - (e_c - c);
}
static double phi_size_mb(const phi *p) {

  if (!p)
    return 0.0;

  double total_mb = (double)sizeof(phi) / (1024.0 * 1024.0);
  uint32_t n = p->n_h;
  for (uint32_t i = 0; i < n; i++) {

    total_mb += c_arr_size_mb(&p->phi_pos[i]);
    total_mb += c_arr_size_mb(&p->phi_inv_pos[i]);
    total_mb += c_arr_size_mb(&p->phi_supp[i]);
    total_mb += c_arr_size_mb(&p->phi_inv_supp[i]);
    total_mb += c_arr_size_mb(&p->phi_l_supp[i]);
  }

  return total_mb;
}
#endif
