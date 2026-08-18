#include "pbwt_col.h"
#include "c_arr.h"
#include "kvec.h"
#include "utils.h"
#include <stdio.h>

pbwt_col build_col(const uint8_t *c, const uint32_t *pa, const uint32_t *da,
                   uint32_t n_h) {
  uint32_t count0 = 0;
  uint32_t count0tmp = 0;
  uint32_t count1 = 0;
  uint32_t u = 0, v = 0;

  int_vec p;
  kv_init(p);
  int_vec uv;
  kv_init(uv);
  int_vec t;
  kv_init(t);
  int_vec b_pa;
  kv_init(b_pa);
  int_vec e_pa;
  kv_init(e_pa);
  int_vec b_da;
  kv_init(b_da);

  bool start = true;

  for (size_t i = 0; i < n_h; ++i) {
    if (i == 0 && c[pa[i]] == 1)
      start = false;
    if (c[i] == 0)
      count0++;
  }

  uint32_t lcs = 0;

  uint32_t p_tmp = 0;
  uint32_t tmp_beg = 0;
  uint32_t tmp_lcp = 0;
  uint32_t tmp_thr = 0;

  bool begrun = true;
  bool pushz = false;
  bool pusho = false;

  if (start)
    pusho = true;
  else
    pushz = true;

  for (size_t i = 0; i < n_h; i++) {
    if (begrun) {
      tmp_beg = pa[i];
      tmp_lcp = da[i];
      u = count0tmp;
      v = count1;
      begrun = false;
    }
    if (c[pa[i]] == 1)
      count1++;
    else
      count0tmp++;

    if (i == 0 || c[pa[i]] != c[pa[i - 1]]) {
      tmp_thr = i;
      lcs = da[i];
    }
    if (da[i] < lcs) {
      tmp_thr = i;
      lcs = da[i];
    }

    if (i == 0 || c[pa[i]] != c[pa[i - 1]]) {
      p_tmp = i;
    }

    if (i == n_h - 1 || c[pa[i]] != c[pa[i + 1]]) {

      if (pusho) {
        kv_push(uint32_t, p, p_tmp);
        kv_push(uint32_t, uv, v);
        bool tmp = pusho;
        pusho = pushz;
        pushz = tmp;
      } else {
        kv_push(uint32_t, p, p_tmp);
        kv_push(uint32_t, uv, u);
        bool tmp = pusho;
        pusho = pushz;
        pushz = tmp;
      }

      if (i + 1 != n_h && da[i + 1] < da[tmp_thr])
        kv_push(uint32_t, t, i + 1);
      else
        kv_push(uint32_t, t, tmp_thr);

      kv_push(uint32_t, b_pa, tmp_beg);
      kv_push(uint32_t, e_pa, pa[i]);
      kv_push(uint32_t, b_da, tmp_lcp);

      begrun = true;
    }
  }
  pbwt_col col;
  col.zero = start;
  col.n_zeros = count0;

  c_arr *tm;

  tm = c_arr_from_kvec(p);
  col.p = *tm;
  free(tm);
  kv_destroy(p);
  tm = c_arr_from_kvec(uv);
  col.uv = *tm;
  free(tm);
  kv_destroy(uv);
  tm = c_arr_from_kvec(t);
  col.t = *tm;
  free(tm);
  kv_destroy(t);
  tm = c_arr_from_kvec(b_pa);
  col.b_pa = *tm;
  free(tm);
  kv_destroy(b_pa);
  tm = c_arr_from_kvec(e_pa);
  col.e_pa = *tm;
  free(tm);
  kv_destroy(e_pa);
  tm = c_arr_from_kvec(b_da);
  col.b_da = *tm;
  free(tm);
  kv_destroy(b_da);

  return col;
}
void print_pbwt_col(const pbwt_col *col, uint32_t n_h) {
  if (!col)
    return;

  printf("Having zero = %d, and zero_c = %d\n", col->zero, col->n_zeros);

  printf("p:\n");
  c_arr_print(&col->p);

  printf("uv:\n");
  c_arr_print(&col->uv);

  printf("t:\n");
  c_arr_print(&col->t);

  printf("b_pa:\n");
  c_arr_print(&col->b_pa);

  printf("e_pa:\n");
  c_arr_print(&col->e_pa);

  printf("b_da:\n");
  c_arr_print(&col->b_da);
}
void free_pbwt_col(pbwt_col *col) {
  if (!col)
    return;

  c_arr_free(&col->p);
  c_arr_free(&col->uv);
  c_arr_free(&col->t);
  c_arr_free(&col->b_pa);
  c_arr_free(&col->b_da);
  c_arr_free(&col->e_pa);
}
int pbwt_col_serialize(FILE *fp, const pbwt_col *col) {
  if (!fp || !col)
    return -1;

  c_arr_serialize(&col->p, fp);
  c_arr_serialize(&col->uv, fp);
  c_arr_serialize(&col->t, fp);
  c_arr_serialize(&col->b_pa, fp);
  c_arr_serialize(&col->e_pa, fp);
  c_arr_serialize(&col->b_da, fp);

  fwrite(&col->zero, sizeof(bool), 1, fp);

  fwrite(&col->n_zeros, sizeof(uint32_t), 1, fp);

  return 0;
}

int pbwt_col_deserialize(FILE *fp, pbwt_col *col) {
  if (!fp || !col)
    return -1;

  c_arr *tmp;

  if (!(tmp = c_arr_deserialize(fp)))
    return -1;
  col->p = *tmp;
  free(tmp);

  if (!(tmp = c_arr_deserialize(fp)))
    return -1;
  col->uv = *tmp;
  free(tmp);

  if (!(tmp = c_arr_deserialize(fp)))
    return -1;
  col->t = *tmp;
  free(tmp);

  if (!(tmp = c_arr_deserialize(fp)))
    return -1;
  col->b_pa = *tmp;
  free(tmp);

  if (!(tmp = c_arr_deserialize(fp)))
    return -1;
  col->e_pa = *tmp;
  free(tmp);

  if (!(tmp = c_arr_deserialize(fp)))
    return -1;
  col->b_da = *tmp;
  free(tmp);

  if (fread(&col->zero, sizeof(bool), 1, fp) != 1)
    return -1;

  if (fread(&col->n_zeros, sizeof(uint32_t), 1, fp) != 1)
    return -1;

  return 0;
}
int pbwt_col_serialize_gz(gzFile fp, const pbwt_col *col) {
  if (!fp || !col)
    return -1;

  c_arr_serialize_gz(fp, &col->p);
  c_arr_serialize_gz(fp, &col->uv);
  c_arr_serialize_gz(fp, &col->t);
  c_arr_serialize_gz(fp, &col->b_pa);
  c_arr_serialize_gz(fp, &col->e_pa);
  c_arr_serialize_gz(fp, &col->b_da);

  gzwrite(fp, &col->zero, sizeof(bool));
  gzwrite(fp, &col->n_zeros, sizeof(uint32_t));

  return 0;
}

int pbwt_col_deserialize_gz(gzFile fp, pbwt_col *col) {
  if (!fp || !col)
    return -1;

  if (c_arr_deserialize_inplace_gz(fp, &col->p) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &col->uv) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &col->t) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &col->b_pa) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &col->e_pa) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &col->b_da) != 0)
    return -1;

  if (gzread(fp, &col->zero, sizeof(bool)) != sizeof(bool))
    return -1;
  if (gzread(fp, &col->n_zeros, sizeof(uint32_t)) != sizeof(uint32_t))
    return -1;

  return 0;
}
