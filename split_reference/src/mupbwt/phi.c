#include "phi.h"
#include "c_arr.h"
#include "utils.h"

static uint32_t *flatten_offsets(const int_vec *arrs, uint32_t n_h) {
  uint32_t *off = (uint32_t *)malloc((n_h + 1) * sizeof(uint32_t));
  off[0] = 0;
  for (uint32_t i = 0; i < n_h; i++)
    off[i + 1] = off[i] + (uint32_t)arrs[i].n;
  return off;
}

static c_arr flatten_field(const int_vec *arrs, uint32_t n_h,
                           const uint32_t *off) {
  uint32_t total = off[n_h];
  if (total == 0)
    return (c_arr){.data = NULL, .n = 0, .bits = 0};

  uint32_t max_val = 0;
  for (uint32_t i = 0; i < n_h; i++)
    for (size_t j = 0; j < arrs[i].n; j++)
      if (arrs[i].a[j] > max_val)
        max_val = arrs[i].a[j];

  c_arr *tmp = c_arr_create(total, max_val);
  size_t k = 0;
  for (uint32_t i = 0; i < n_h; i++)
    for (size_t j = 0; j < arrs[i].n; j++)
      c_arr_set(tmp, k++, arrs[i].a[j]);

  c_arr result = *tmp;
  free(tmp);
  return result;
}

void build_phi(phi *phi, int_vec *supp_b, int_vec *supp_e, int_vec *supp_pa_b,
               int_vec *supp_pa_e, int_vec *supp_da_b, uint32_t n_h,
               uint32_t n_w) {

  phi->n_h = n_h;
  phi->n_w = n_w;

  phi->phi_pos_off = flatten_offsets(supp_b, n_h);
  phi->phi_pos = flatten_field(supp_b, n_h, phi->phi_pos_off);

  phi->phi_inv_pos_off = flatten_offsets(supp_e, n_h);
  phi->phi_inv_pos = flatten_field(supp_e, n_h, phi->phi_inv_pos_off);

  phi->phi_supp_off = flatten_offsets(supp_pa_b, n_h);
  phi->phi_supp = flatten_field(supp_pa_b, n_h, phi->phi_supp_off);

  phi->phi_inv_supp_off = flatten_offsets(supp_pa_e, n_h);
  phi->phi_inv_supp = flatten_field(supp_pa_e, n_h, phi->phi_inv_supp_off);

  phi->phi_l_supp_off = flatten_offsets(supp_da_b, n_h);
  phi->phi_l_supp = flatten_field(supp_da_b, n_h, phi->phi_l_supp_off);
}

static void print_field(const c_arr *data, const uint32_t *off, uint32_t n_h) {
  for (uint32_t i = 0; i < n_h; i++) {
    uint32_t base = off[i];
    uint32_t len = off[i + 1] - base;
    for (uint32_t j = 0; j < len; j++)
      printf("%u ", c_arr_get(data, base + j));
    printf("\n");
  }
}

void phi_print(const phi *p) {
  if (!p)
    return;
  printf("phi_pos\n");
  print_field(&p->phi_pos, p->phi_pos_off, p->n_h);
  printf("phi_inv_pos\n");
  print_field(&p->phi_inv_pos, p->phi_inv_pos_off, p->n_h);
  printf("phi_supp\n");
  print_field(&p->phi_supp, p->phi_supp_off, p->n_h);
  printf("phi_inv_supp\n");
  print_field(&p->phi_inv_supp, p->phi_inv_supp_off, p->n_h);
  printf("phi_l_supp\n");
  print_field(&p->phi_l_supp, p->phi_l_supp_off, p->n_h);
}

void phi_free(phi *p) {
  if (!p)
    return;

  c_arr_free(&p->phi_pos);
  c_arr_free(&p->phi_inv_pos);
  c_arr_free(&p->phi_supp);
  c_arr_free(&p->phi_inv_supp);
  c_arr_free(&p->phi_l_supp);

  free(p->phi_pos_off);
  p->phi_pos_off = NULL;
  free(p->phi_inv_pos_off);
  p->phi_inv_pos_off = NULL;
  free(p->phi_supp_off);
  p->phi_supp_off = NULL;
  free(p->phi_inv_supp_off);
  p->phi_inv_supp_off = NULL;
  free(p->phi_l_supp_off);
  p->phi_l_supp_off = NULL;

  p->n_h = 0;
  p->n_w = 0;
}

static void write_offsets(const uint32_t *off, uint32_t n_h, FILE *fp) {
  fwrite(off, sizeof(uint32_t), n_h + 1, fp);
}

static void read_offsets(uint32_t **off, uint32_t n_h, FILE *fp) {
  *off = (uint32_t *)malloc((n_h + 1) * sizeof(uint32_t));
  fread(*off, sizeof(uint32_t), n_h + 1, fp);
}

int phi_serialize(const phi *p, FILE *fp) {
  if (!fp || !p)
    return -1;

  fwrite(&p->n_h, sizeof(uint32_t), 1, fp);
  fwrite(&p->n_w, sizeof(uint32_t), 1, fp);

  write_offsets(p->phi_pos_off, p->n_h, fp);
  c_arr_serialize(&p->phi_pos, fp);
  write_offsets(p->phi_inv_pos_off, p->n_h, fp);
  c_arr_serialize(&p->phi_inv_pos, fp);
  write_offsets(p->phi_supp_off, p->n_h, fp);
  c_arr_serialize(&p->phi_supp, fp);
  write_offsets(p->phi_inv_supp_off, p->n_h, fp);
  c_arr_serialize(&p->phi_inv_supp, fp);
  write_offsets(p->phi_l_supp_off, p->n_h, fp);
  c_arr_serialize(&p->phi_l_supp, fp);

  return 0;
}

int phi_deserialize(phi *p, FILE *fp) {
  if (!fp || !p)
    return -1;

  fread(&p->n_h, sizeof(uint32_t), 1, fp);
  fread(&p->n_w, sizeof(uint32_t), 1, fp);

  c_arr *tmp;

  read_offsets(&p->phi_pos_off, p->n_h, fp);
  tmp = c_arr_deserialize(fp);
  p->phi_pos = *tmp;
  free(tmp);

  read_offsets(&p->phi_inv_pos_off, p->n_h, fp);
  tmp = c_arr_deserialize(fp);
  p->phi_inv_pos = *tmp;
  free(tmp);

  read_offsets(&p->phi_supp_off, p->n_h, fp);
  tmp = c_arr_deserialize(fp);
  p->phi_supp = *tmp;
  free(tmp);

  read_offsets(&p->phi_inv_supp_off, p->n_h, fp);
  tmp = c_arr_deserialize(fp);
  p->phi_inv_supp = *tmp;
  free(tmp);

  read_offsets(&p->phi_l_supp_off, p->n_h, fp);
  tmp = c_arr_deserialize(fp);
  p->phi_l_supp = *tmp;
  free(tmp);

  return 0;
}

static void write_offsets_gz(const uint32_t *off, uint32_t n_h, gzFile fp) {
  gzwrite(fp, off, (unsigned int)((n_h + 1) * sizeof(uint32_t)));
}

static int read_offsets_gz(uint32_t **off, uint32_t n_h, gzFile fp) {
  *off = (uint32_t *)malloc((n_h + 1) * sizeof(uint32_t));
  unsigned int want = (unsigned int)((n_h + 1) * sizeof(uint32_t));
  if (gzread(fp, *off, want) != (int)want)
    return -1;
  return 0;
}

int phi_serialize_gz(const phi *p, gzFile fp) {
  if (!fp || !p)
    return -1;

  gzwrite(fp, &p->n_h, sizeof(uint32_t));
  gzwrite(fp, &p->n_w, sizeof(uint32_t));

  write_offsets_gz(p->phi_pos_off, p->n_h, fp);
  c_arr_serialize_gz(fp, &p->phi_pos);
  write_offsets_gz(p->phi_inv_pos_off, p->n_h, fp);
  c_arr_serialize_gz(fp, &p->phi_inv_pos);
  write_offsets_gz(p->phi_supp_off, p->n_h, fp);
  c_arr_serialize_gz(fp, &p->phi_supp);
  write_offsets_gz(p->phi_inv_supp_off, p->n_h, fp);
  c_arr_serialize_gz(fp, &p->phi_inv_supp);
  write_offsets_gz(p->phi_l_supp_off, p->n_h, fp);
  c_arr_serialize_gz(fp, &p->phi_l_supp);

  return 0;
}

int phi_deserialize_gz(phi *p, gzFile fp) {
  if (!fp || !p)
    return -1;

  if (gzread(fp, &p->n_h, sizeof(uint32_t)) != sizeof(uint32_t))
    return -1;
  if (gzread(fp, &p->n_w, sizeof(uint32_t)) != sizeof(uint32_t))
    return -1;

  if (read_offsets_gz(&p->phi_pos_off, p->n_h, fp) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &p->phi_pos) != 0)
    return -1;

  if (read_offsets_gz(&p->phi_inv_pos_off, p->n_h, fp) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &p->phi_inv_pos) != 0)
    return -1;

  if (read_offsets_gz(&p->phi_supp_off, p->n_h, fp) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &p->phi_supp) != 0)
    return -1;

  if (read_offsets_gz(&p->phi_inv_supp_off, p->n_h, fp) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &p->phi_inv_supp) != 0)
    return -1;

  if (read_offsets_gz(&p->phi_l_supp_off, p->n_h, fp) != 0)
    return -1;
  if (c_arr_deserialize_inplace_gz(fp, &p->phi_l_supp) != 0)
    return -1;

  return 0;
}
