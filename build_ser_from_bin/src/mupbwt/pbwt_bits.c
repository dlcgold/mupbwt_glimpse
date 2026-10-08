#include "pbwt.h"
#include "pbwt_col.h"
#include "kvec.h"
#include "c_arr.h"
#include <string.h>

void pbwt_build_af_process_site_from_bits(pbwt *pbwt, pbwt_build_state *st,
                                          const uint8_t *col_bits) {
  pbwt->n_sites++;

  memcpy(st->c_col, col_bits, (size_t)pbwt->n_haps * sizeof(uint8_t));

  pbwt_col col = build_col(st->c_col, st->pa, st->da, pbwt->n_haps);
  for (size_t r = 0; r < col.p.n; r++) {
    kv_push(uint32_t, st->supp_b[c_arr_get(&col.b_pa, r)], st->c);
    kv_push(uint32_t, st->supp_e[c_arr_get(&col.e_pa, r)], st->c);
    if (r == 0) {
      kv_push(uint32_t, st->supp_pa_b[c_arr_get(&col.b_pa, r)], pbwt->n_haps);
      kv_push(uint32_t, st->supp_da_b[c_arr_get(&col.b_pa, r)], 0);
    } else {
      kv_push(uint32_t, st->supp_pa_b[c_arr_get(&col.b_pa, r)],
              st->pa[c_arr_get(&col.p, r) - 1]);
      kv_push(uint32_t, st->supp_da_b[c_arr_get(&col.b_pa, r)],
              st->da[c_arr_get(&col.p, r)]);
    }
    if (r == col.p.n - 1) {
      kv_push(uint32_t, st->supp_pa_e[c_arr_get(&col.e_pa, r)], pbwt->n_haps);
    } else {
      kv_push(uint32_t, st->supp_pa_e[c_arr_get(&col.e_pa, r)],
              st->pa[c_arr_get(&col.p, r + 1)]);
    }
  }

  memcpy(st->l_pa, st->pa, pbwt->n_haps * sizeof(uint32_t));
  memcpy(st->l_da, st->da, pbwt->n_haps * sizeof(uint32_t));

  KV_PUSH_COL_VEC(pbwt->cols, col);

  pbwt_update(st->c_col, &st->pa, &st->da, pbwt->n_haps);
  st->c++;
}
