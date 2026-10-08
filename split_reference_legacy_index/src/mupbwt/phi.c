#include "phi.h"
#include "c_arr.h"
#include "utils.h"

void build_phi(phi *phi, int_vec *supp_b, int_vec *supp_e, int_vec *supp_pa_b,
               int_vec *supp_pa_e, int_vec *supp_da_b, uint32_t n_h,
               uint32_t n_w) {

  phi->n_h = n_h;
  phi->n_w = n_w;

  phi->phi_pos = (c_arr *)malloc(n_h * sizeof(c_arr));
  phi->phi_inv_pos = (c_arr *)malloc(n_h * sizeof(c_arr));
  phi->phi_supp = (c_arr *)malloc(n_h * sizeof(c_arr));
  phi->phi_inv_supp = (c_arr *)malloc(n_h * sizeof(c_arr));
  phi->phi_l_supp = (c_arr *)malloc(n_h * sizeof(c_arr));

  for (size_t i = 0; i < n_h; i++) {

    if (supp_b[i].n > 0) {
      c_arr *s_b = c_arr_from_kvec(supp_b[i]);
      phi->phi_pos[i] = *s_b;
      free(s_b);
    } else {
      phi->phi_pos[i] = (c_arr){.data = NULL, .n = 0, .bits = 0};
    }

    if (supp_e[i].n > 0) {
      c_arr *s_e = c_arr_from_kvec(supp_e[i]);
      phi->phi_inv_pos[i] = *s_e;
      free(s_e);
    } else {
      phi->phi_inv_pos[i] = (c_arr){.data = NULL, .n = 0, .bits = 0};
    }

    if (supp_pa_b[i].n > 0) {
      c_arr *s_pb = c_arr_from_kvec(supp_pa_b[i]);
      phi->phi_supp[i] = *s_pb;
      free(s_pb);
    } else {
      phi->phi_supp[i] = (c_arr){.data = NULL, .n = 0, .bits = 0};
    }

    if (supp_pa_e[i].n > 0) {
      c_arr *s_pe = c_arr_from_kvec(supp_pa_e[i]);
      phi->phi_inv_supp[i] = *s_pe;
      free(s_pe);
    } else {
      phi->phi_inv_supp[i] = (c_arr){.data = NULL, .n = 0, .bits = 0};
    }

    if (supp_da_b[i].n > 0) {
      c_arr *l_b = c_arr_from_kvec(supp_da_b[i]);
      phi->phi_l_supp[i] = *l_b;
      free(l_b);
    } else {
      phi->phi_l_supp[i] = (c_arr){.data = NULL, .n = 0, .bits = 0};
    }
  }
}

void phi_print(const phi *p) {
  if (!p)
    return;
  for (size_t i = 0; i < p->n_h; i++) {
    printf("phi stuff for row %zu\n", i);
    if (p->phi_pos[i].n > 0)
      c_arr_print(&p->phi_pos[i]);
    else
      printf("\n");

    if (p->phi_inv_pos[i].n > 0)
      c_arr_print(&p->phi_inv_pos[i]);
    else
      printf("\n");

    if (p->phi_supp[i].n > 0)
      c_arr_print(&p->phi_supp[i]);
    else
      printf("\n");

    if (p->phi_inv_supp[i].n > 0)
      c_arr_print(&p->phi_inv_supp[i]);
    else
      printf("\n");

    if (p->phi_l_supp[i].n > 0)
      c_arr_print(&p->phi_l_supp[i]);
    else
      printf("\n");
    printf("\n");
  }
}

void phi_free(phi *p) {
  if (!p)
    return;

  for (uint32_t i = 0; i < p->n_h; i++) {
    c_arr_free(&p->phi_pos[i]);
    c_arr_free(&p->phi_inv_pos[i]);
    c_arr_free(&p->phi_supp[i]);
    c_arr_free(&p->phi_inv_supp[i]);
    c_arr_free(&p->phi_l_supp[i]);
  }

  free(p->phi_pos);
  p->phi_pos = NULL;
  free(p->phi_inv_pos);
  p->phi_inv_pos = NULL;
  free(p->phi_supp);
  p->phi_supp = NULL;
  free(p->phi_inv_supp);
  p->phi_inv_supp = NULL;
  free(p->phi_l_supp);
  p->phi_l_supp = NULL;

  p->n_h = 0;
  p->n_w = 0;
}

int phi_serialize(const phi *p, FILE *fp) {
  if (!fp || !p)
    return -1;

  fwrite(&p->n_h, sizeof(uint32_t), 1, fp);
  fwrite(&p->n_w, sizeof(uint32_t), 1, fp);

  for (uint32_t i = 0; i < p->n_h; i++) {
    c_arr_serialize(&p->phi_pos[i], fp);
    c_arr_serialize(&p->phi_inv_pos[i], fp);
    c_arr_serialize(&p->phi_supp[i], fp);
    c_arr_serialize(&p->phi_inv_supp[i], fp);
    c_arr_serialize(&p->phi_l_supp[i], fp);
  }

  return 0;
}

int phi_deserialize(phi *p, FILE *fp) {
  if (!fp || !p)
    return -1;

  fread(&p->n_h, sizeof(uint32_t), 1, fp);
  fread(&p->n_w, sizeof(uint32_t), 1, fp);

  p->phi_pos = (c_arr *)malloc(p->n_h * sizeof(c_arr));
  p->phi_inv_pos = (c_arr *)malloc(p->n_h * sizeof(c_arr));
  p->phi_supp = (c_arr *)malloc(p->n_h * sizeof(c_arr));
  p->phi_inv_supp = (c_arr *)malloc(p->n_h * sizeof(c_arr));
  p->phi_l_supp = (c_arr *)malloc(p->n_h * sizeof(c_arr));

  for (uint32_t i = 0; i < p->n_h; i++) {
    c_arr *tmp;

    tmp = c_arr_deserialize(fp);
    p->phi_pos[i] = *tmp;
    free(tmp);

    tmp = c_arr_deserialize(fp);
    p->phi_inv_pos[i] = *tmp;
    free(tmp);

    tmp = c_arr_deserialize(fp);
    p->phi_supp[i] = *tmp;
    free(tmp);

    tmp = c_arr_deserialize(fp);
    p->phi_inv_supp[i] = *tmp;
    free(tmp);

    tmp = c_arr_deserialize(fp);
    p->phi_l_supp[i] = *tmp;
    free(tmp);
  }

  return 0;
}

int phi_serialize_gz(const phi *p, gzFile fp) {
  if (!fp || !p)
    return -1;

  gzwrite(fp, &p->n_h, sizeof(uint32_t));
  gzwrite(fp, &p->n_w, sizeof(uint32_t));

  for (uint32_t i = 0; i < p->n_h; i++) {
    c_arr_serialize_gz(fp, &p->phi_pos[i]);
    c_arr_serialize_gz(fp, &p->phi_inv_pos[i]);
    c_arr_serialize_gz(fp, &p->phi_supp[i]);
    c_arr_serialize_gz(fp, &p->phi_inv_supp[i]);
    c_arr_serialize_gz(fp, &p->phi_l_supp[i]);
  }

  return 0;
}
int phi_deserialize_gz(phi *p, gzFile fp) {
  if (!fp || !p)
    return -1;

  if (gzread(fp, &p->n_h, sizeof(uint32_t)) != sizeof(uint32_t))
    return -1;
  if (gzread(fp, &p->n_w, sizeof(uint32_t)) != sizeof(uint32_t))
    return -1;

  p->phi_pos = (c_arr *)malloc(p->n_h * sizeof(c_arr));
  p->phi_inv_pos = (c_arr *)malloc(p->n_h * sizeof(c_arr));
  p->phi_supp = (c_arr *)malloc(p->n_h * sizeof(c_arr));
  p->phi_inv_supp = (c_arr *)malloc(p->n_h * sizeof(c_arr));
  p->phi_l_supp = (c_arr *)malloc(p->n_h * sizeof(c_arr));

  if (!p->phi_pos || !p->phi_inv_pos || !p->phi_supp || !p->phi_inv_supp ||
      !p->phi_l_supp) {
    return -1;
  }

  for (uint32_t i = 0; i < p->n_h; i++) {
    if (c_arr_deserialize_inplace_gz(fp, &p->phi_pos[i]) != 0)
      return -1;
    if (c_arr_deserialize_inplace_gz(fp, &p->phi_inv_pos[i]) != 0)
      return -1;
    if (c_arr_deserialize_inplace_gz(fp, &p->phi_supp[i]) != 0)
      return -1;
    if (c_arr_deserialize_inplace_gz(fp, &p->phi_inv_supp[i]) != 0)
      return -1;
    if (c_arr_deserialize_inplace_gz(fp, &p->phi_l_supp[i]) != 0)
      return -1;
  }

  return 0;
}
