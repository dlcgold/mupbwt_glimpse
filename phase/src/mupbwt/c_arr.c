#include "c_arr.h"

static uint8_t calc_bits(uint32_t max_val) {
  return max_val == 0 ? 1 : (uint8_t)(floor(log2(max_val)) + 1);
}

c_arr *c_arr_create(size_t n, uint32_t max_val) {
  c_arr *arr = (c_arr *)malloc(sizeof(c_arr));
  if (!arr)
    return NULL;
  arr->n = n;
  arr->bits = calc_bits(max_val);
  size_t total_bits = n * arr->bits;
  size_t uint32_count = (total_bits + 31) / 32;
  arr->data = (uint32_t *)calloc(uint32_count, sizeof(uint32_t));
  if (!arr->data) {
    free(arr);
    return NULL;
  }
  return arr;
}

c_arr *c_arr_from_array(const uint32_t *src, size_t n) {
  if (n == 0)
    return NULL;
  uint32_t max_val = 0;
  for (size_t i = 0; i < n; ++i)
    if (src[i] > max_val)
      max_val = src[i];
  c_arr *arr = c_arr_create(n, max_val);
  if (!arr)
    return NULL;
  for (size_t i = 0; i < n; ++i)
    c_arr_set(arr, i, src[i]);
  return arr;
}

c_arr *c_arr_from_kvec(const int_vec kv) {
  return c_arr_from_array(kv.a, kv.n);
}

c_arr *c_arr_from_array_max(const uint32_t *src, size_t n, uint32_t max_val) {
  if (n == 0)
    return NULL;
  c_arr *arr = c_arr_create(n, max_val);
  if (!arr)
    return NULL;
  for (size_t i = 0; i < n; ++i)
    c_arr_set(arr, i, src[i]);
  return arr;
}

c_arr *c_arr_from_kvec_max(const int_vec kv, uint32_t max_val) {
  return c_arr_from_array_max(kv.a, kv.n, max_val);
}

void c_arr_free(c_arr *arr) {
  if (!arr)
    return;
  if (arr->data)
    free(arr->data);
  arr->data = NULL;
  arr->n = 0;
  arr->bits = 0;
  // free(arr);
}

void c_arr_set(c_arr *arr, size_t i, uint32_t val) {
  if (i >= arr->n)
    return;
  size_t bit_pos = i * arr->bits;
  size_t idx = bit_pos / 32;
  size_t offset = bit_pos % 32;
  arr->data[idx] &= ~(((1u << arr->bits) - 1) << offset);
  arr->data[idx] |= (val & ((1u << arr->bits) - 1)) << offset;
  if (offset + arr->bits > 32) {
    int bits_left = (offset + arr->bits) - 32;
    arr->data[idx + 1] &= ~((1u << bits_left) - 1);
    arr->data[idx + 1] |= val >> (arr->bits - bits_left);
  }
}

bool c_arr_equal(const c_arr *a, const c_arr *b) {
  if (a->n != b->n) {
    return false;
  }

  for (size_t i = 0; i < a->n; i++) {
    if (c_arr_get(a, i) != c_arr_get(b, i))
      return false;
  }
  return true;
}
void c_arr_print(const c_arr *a) {
  for (size_t i = 0; i < a->n; i++) {
    printf("%d ", c_arr_get(a, i));
  }
  printf("\n");
}

uint32_t c_arr_rank(const c_arr *a, uint32_t v) {
  if (!a || a->n == 0)
    return 0;

  size_t pos = 0;
  size_t n = a->n;
  size_t step = 1ULL << (31 - __builtin_clz((uint32_t)n));

  if (c_arr_get(a, 0) >= v)
    return 0;

  for (; step > 0; step >>= 1) {
    size_t next_pos = pos + step;
    if (next_pos < n) {
      if (c_arr_get(a, next_pos) < v) {
        pos = next_pos;
      }
    }
  }

  return (uint32_t)(pos + 1);
}

/* int c_arr_serialize(const c_arr *arr, FILE *fp) { */
/*   if (!arr || !fp) */
/*     return -1; */
/**/
/*   uint64_t n = arr->n; */
/*   uint8_t bits = arr->bits; */
/*   uint64_t total_bits = n * bits; */
/*   uint64_t data_len = (total_bits + 31) / 32; */
/**/
/*   if (fwrite(&n, sizeof(uint64_t), 1, fp) != 1) */
/*     return -1; */
/*   if (fwrite(&bits, sizeof(uint8_t), 1, fp) != 1) */
/*     return -1; */
/*   if (fwrite(&data_len, sizeof(uint64_t), 1, fp) != 1) */
/*     return -1; */
/*   if (fwrite(arr->data, sizeof(uint32_t), data_len, fp) != data_len) */
/*     return -1; */
/**/
/*   return 0; */
/* } */
/**/
/* c_arr *c_arr_deserialize(FILE *fp) { */
/*   if (!fp) */
/*     return NULL; */
/**/
/*   uint64_t n; */
/*   uint8_t bits; */
/*   uint64_t data_len; */
/**/
/*   if (fread(&n, sizeof(uint64_t), 1, fp) != 1) */
/*     return NULL; */
/*   if (fread(&bits, sizeof(uint8_t), 1, fp) != 1) */
/*     return NULL; */
/*   if (fread(&data_len, sizeof(uint64_t), 1, fp) != 1) */
/*     return NULL; */
/**/
/*   c_arr *arr = (c_arr *)malloc(sizeof(c_arr)); */
/*   if (!arr) */
/*     return NULL; */
/**/
/*   arr->n = (size_t)n; */
/*   arr->bits = bits; */
/*   arr->data = (uint32_t *)malloc(data_len * sizeof(uint32_t)); */
/*   if (!arr->data) { */
/*     free(arr); */
/*     return NULL; */
/*   } */
/**/
/*   if (fread(arr->data, sizeof(uint32_t), data_len, fp) != data_len) { */
/*     free(arr->data); */
/*     free(arr); */
/*     return NULL; */
/*   } */
/**/
/*   return arr; */
/* } */

int c_arr_serialize(const c_arr *arr, FILE *fp) {
  if (!arr || !fp)
    return -1;

  uint64_t n = arr->n;
  uint8_t bits = arr->bits;
  uint64_t total_bits = n * bits;

  uint64_t bytes_to_write = (total_bits + 7) / 8;

  if (fwrite(&n, sizeof(uint64_t), 1, fp) != 1)
    return -1;
  if (fwrite(&bits, sizeof(uint8_t), 1, fp) != 1)
    return -1;

  if (fwrite(&bytes_to_write, sizeof(uint64_t), 1, fp) != 1)
    return -1;

  if (fwrite(arr->data, 1, bytes_to_write, fp) != bytes_to_write)
    return -1;

  return 0;
}
c_arr *c_arr_deserialize(FILE *fp) {
  if (!fp)
    return NULL;

  uint64_t n;
  uint8_t bits;
  uint64_t bytes_len;

  if (fread(&n, sizeof(uint64_t), 1, fp) != 1)
    return NULL;
  if (fread(&bits, sizeof(uint8_t), 1, fp) != 1)
    return NULL;
  if (fread(&bytes_len, sizeof(uint64_t), 1, fp) != 1)
    return NULL;

  c_arr *arr = (c_arr *)malloc(sizeof(c_arr));
  if (!arr)
    return NULL;

  arr->n = (size_t)n;
  arr->bits = bits;

  uint64_t n_words = (n * bits + 31) / 32;

  arr->data = (uint32_t *)calloc(n_words + 1, sizeof(uint32_t));
  if (!arr->data) {
    free(arr);
    return NULL;
  }

  if (fread(arr->data, 1, bytes_len, fp) != bytes_len) {
    free(arr->data);
    free(arr);
    return NULL;
  }

  return arr;
}

int c_arr_serialize_gz(gzFile fp, const c_arr *arr) {
  uint64_t n = arr->n;
  uint8_t bits = arr->bits;
  uint64_t bytes_to_write = ((n * bits) + 7) / 8;

  gzwrite(fp, &n, sizeof(uint64_t));
  gzwrite(fp, &bits, sizeof(uint8_t));
  gzwrite(fp, &bytes_to_write, sizeof(uint64_t));
  if (bytes_to_write > 0) {
    gzwrite(fp, arr->data, (unsigned int)bytes_to_write);
  }
  return 0;
}

int c_arr_deserialize_inplace_gz(gzFile fp, c_arr *arr) {
  uint64_t n, bytes_len;
  uint8_t bits;

  if (gzread(fp, &n, sizeof(uint64_t)) != sizeof(uint64_t))
    return -1;
  if (gzread(fp, &bits, sizeof(uint8_t)) != sizeof(uint8_t))
    return -1;
  if (gzread(fp, &bytes_len, sizeof(uint64_t)) != sizeof(uint64_t))
    return -1;

  arr->n = (size_t)n;
  arr->bits = bits;

  uint64_t n_words = (n * bits + 31) / 32;
  arr->data = (uint32_t *)calloc(n_words + 1, sizeof(uint32_t));

  if (bytes_len > 0) {
    if (gzread(fp, arr->data, (unsigned int)bytes_len) != (int)bytes_len)
      return -1;
  }
  return 0;
}
