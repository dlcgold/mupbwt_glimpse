#ifndef PHI_H
#define PHI_H

#include "c_arr.h"
#include "kvec.h"
#include "utils.h"

typedef struct {
  c_arr pos;
  uint32_t *pos_off;
  c_arr inv_pos;
  uint32_t *inv_pos_off;
  c_arr supp;
  uint32_t *supp_off;
  c_arr inv_supp;
  uint32_t *inv_supp_off;
  c_arr l_supp;
  uint32_t *l_supp_off;
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

static inline __attribute__((always_inline)) uint32_t phi_f(const phi *p,
                                                            uint32_t pa_v,
                                                            uint32_t c) {
  uint32_t pos_base = p->pos_off[pa_v];
  uint32_t pos_len = p->pos_off[pa_v + 1] - pos_base;
  uint32_t supp_base = p->supp_off[pa_v];
  uint32_t supp_len = p->supp_off[pa_v + 1] - supp_base;

  uint32_t t_c = c_arr_rank_range(&p->pos, pos_base, pos_len, c);
  t_c -= (t_c == supp_len) & (t_c > 0);
  return c_arr_get(&p->supp, supp_base + t_c);
}

static inline __attribute__((always_inline)) uint32_t phi_inv_f(const phi *p,
                                                                uint32_t pa_v,
                                                                uint32_t c) {
  uint32_t pos_base = p->inv_pos_off[pa_v];
  uint32_t pos_len = p->inv_pos_off[pa_v + 1] - pos_base;
  uint32_t supp_base = p->inv_supp_off[pa_v];
  uint32_t supp_len = p->inv_supp_off[pa_v + 1] - supp_base;

  uint32_t t_c = c_arr_rank_range(&p->inv_pos, pos_base, pos_len, c);
  t_c -= (t_c == supp_len) & (t_c > 0);
  return c_arr_get(&p->inv_supp, supp_base + t_c);
}

static inline __attribute__((always_inline)) uint32_t phi_l(const phi *p,
                                                            uint32_t pa_v,
                                                            uint32_t c) {
  if (c == 0 || phi_f(p, pa_v, c) == p->n_h)
    return 0;

  uint32_t pos_base = p->pos_off[pa_v];
  uint32_t pos_len = p->pos_off[pa_v + 1] - pos_base;
  uint32_t supp_base = p->supp_off[pa_v];
  uint32_t supp_len = p->supp_off[pa_v + 1] - supp_base;
  uint32_t l_supp_base = p->l_supp_off[pa_v];

  uint32_t t_c = c_arr_rank_range(&p->pos, pos_base, pos_len, c);
  t_c -= (t_c == supp_len) & (t_c > 0);

  uint32_t e_c = (t_c < pos_len) ? c_arr_get(&p->pos, pos_base + t_c) : p->n_w - 1;

  return c_arr_get(&p->l_supp, l_supp_base + t_c) - (e_c - c);
}

static double phi_size_mb(const phi *p) {

  if (!p)
    return 0.0;

  double total_mb = (double)sizeof(phi) / (1024.0 * 1024.0);
  total_mb += c_arr_size_mb(&p->pos);
  total_mb += c_arr_size_mb(&p->inv_pos);
  total_mb += c_arr_size_mb(&p->supp);
  total_mb += c_arr_size_mb(&p->inv_supp);
  total_mb += c_arr_size_mb(&p->l_supp);
  total_mb += 5.0 * (p->n_h + 1) * sizeof(uint32_t) / (1024.0 * 1024.0);

  return total_mb;
}
#endif
