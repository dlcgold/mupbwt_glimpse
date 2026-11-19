#include "phi_fix.h"
#include "c_arr.h"
#include <stdlib.h>

int phi_fixup(phi *p) {
  uint32_t n_h = p->n_h, last = p->n_w - 1;
  uint32_t n_fix = 0;

  for (uint32_t h = 0; h < n_h; h++) {
    uint32_t pl = p->phi_pos_off[h + 1] - p->phi_pos_off[h];
    uint32_t sl = p->phi_supp_off[h + 1] - p->phi_supp_off[h];
    if (pl == 0 || sl != pl)
      continue;
    uint32_t e = c_arr_get(&p->phi_pos, p->phi_pos_off[h] + pl - 1);
    uint32_t pred = c_arr_get(&p->phi_supp, p->phi_supp_off[h] + sl - 1);
    if (e < last && pred != n_h)
      n_fix++;
  }
  if (n_fix == 0)
    return 0;

  size_t n_supp = p->phi_supp.n + n_fix, n_l = p->phi_l_supp.n + n_fix;
  uint32_t max_l = last;
  c_arr *supp = c_arr_create(n_supp, n_h);
  c_arr *l_supp = c_arr_create(n_l, max_l);
  uint32_t *supp_off = (uint32_t *)malloc((n_h + 1) * sizeof(uint32_t));
  uint32_t *l_supp_off = (uint32_t *)malloc((n_h + 1) * sizeof(uint32_t));
  if (!supp || !l_supp || !supp_off || !l_supp_off) {
    free(supp_off);
    free(l_supp_off);
    return -1;
  }

  size_t k = 0, kl = 0;
  supp_off[0] = l_supp_off[0] = 0;
  for (uint32_t h = 0; h < n_h; h++) {
    uint32_t pl = p->phi_pos_off[h + 1] - p->phi_pos_off[h];
    uint32_t sb = p->phi_supp_off[h], sl = p->phi_supp_off[h + 1] - sb;
    uint32_t lb = p->phi_l_supp_off[h], ll = p->phi_l_supp_off[h + 1] - lb;
    for (uint32_t j = 0; j < sl; j++)
      c_arr_set(supp, k++, c_arr_get(&p->phi_supp, sb + j));
    for (uint32_t j = 0; j < ll; j++)
      c_arr_set(l_supp, kl++, c_arr_get(&p->phi_l_supp, lb + j));
    if (pl > 0 && sl == pl) {
      uint32_t e = c_arr_get(&p->phi_pos, p->phi_pos_off[h] + pl - 1);
      uint32_t pred = c_arr_get(&p->phi_supp, sb + sl - 1);
      if (e < last && pred != n_h) {
        c_arr_set(supp, k++, pred);
        c_arr_set(l_supp, kl++, last - 1 - e);
      }
    }
    supp_off[h + 1] = (uint32_t)k;
    l_supp_off[h + 1] = (uint32_t)kl;
  }

  c_arr_free(&p->phi_supp);
  c_arr_free(&p->phi_l_supp);
  free(p->phi_supp_off);
  free(p->phi_l_supp_off);
  p->phi_supp = *supp;
  p->phi_l_supp = *l_supp;
  free(supp);
  free(l_supp);
  p->phi_supp_off = supp_off;
  p->phi_l_supp_off = l_supp_off;
  return (int)n_fix;
}
