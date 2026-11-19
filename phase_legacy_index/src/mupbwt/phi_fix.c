#include "phi_fix.h"
#include "c_arr.h"
#include <stdlib.h>

static int append_entry(c_arr *a, uint32_t val) {
  size_t n = a->n;
  uint32_t max_val = val;
  for (size_t j = 0; j < n; j++) {
    uint32_t v = c_arr_get(a, j);
    if (v > max_val)
      max_val = v;
  }
  c_arr *b = c_arr_create(n + 1, max_val);
  if (!b)
    return -1;
  for (size_t j = 0; j < n; j++)
    c_arr_set(b, j, c_arr_get(a, j));
  c_arr_set(b, n, val);
  c_arr_free(a);
  *a = *b;
  free(b);
  return 0;
}

int phi_fixup(phi *p) {
  uint32_t n_h = p->n_h, last = p->n_w - 1;
  uint32_t n_fix = 0;

  for (uint32_t h = 0; h < n_h; h++) {
    size_t pl = p->phi_pos[h].n;
    size_t sl = p->phi_supp[h].n;
    if (pl == 0 || sl != pl)
      continue;
    uint32_t e = c_arr_get(&p->phi_pos[h], pl - 1);
    uint32_t pred = c_arr_get(&p->phi_supp[h], sl - 1);
    if (e < last && pred != n_h) {
      if (append_entry(&p->phi_supp[h], pred) != 0)
        return -1;
      if (append_entry(&p->phi_l_supp[h], last - 1 - e) != 0)
        return -1;
      n_fix++;
    }
  }
  return (int)n_fix;
}
