#include "pbwt.h"
#include "c_arr.h"
#include "kvec.h"
#include "match.h"
#include "pbwt_col.h"
#include "phi.h"
#include "utils.h"
#include <cstdint>
#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

/* void pbwt_build(char *filename, pbwt *pbwt, int threads) { */
/*   omp_set_num_threads(threads); */
/*   htsFile *fp = hts_open(filename, "rb"); */
/**/
/*   fprintf(stderr, "Reading %s\n", filename); */
/**/
/*   if (fp == NULL) { */
/*     fprintf(stderr, "%s NOT FOUND\n", filename); */
/*     exit(-1); */
/*   } */
/**/
/*   bcf_hdr_t *hdr = bcf_hdr_read(fp); */
/*   bcf1_t *rec = bcf_init(); */
/**/
/*   pbwt->n_haps = bcf_hdr_nsamples(hdr) * 2; */
/*   pbwt->n_sites = 0; */
/**/
/*   fprintf(stderr, "# haplotypes = %d\n", pbwt->n_haps); */
/**/
/*   uint32_t *pa = (uint32_t *)malloc(pbwt->n_haps * sizeof(uint32_t)); */
/*   uint32_t *da = (uint32_t *)malloc(pbwt->n_haps * sizeof(uint32_t)); */
/*   uint32_t *l_pa = (uint32_t *)malloc(pbwt->n_haps * sizeof(uint32_t)); */
/*   uint32_t *l_da = (uint32_t *)malloc(pbwt->n_haps * sizeof(uint32_t)); */
/**/
/*   // uint32_t *thr = (uint32_t *)malloc(pbwt.n_haps * sizeof(uint32_t)); */
/*   // uint32_t *thr = NULL; */
/**/
/*   int_vec *supp_b = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec)); */
/*   int_vec *supp_e = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec)); */
/**/
/*   int_vec *supp_pa_b = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec)); */
/*   int_vec *supp_da_b = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec)); */
/*   int_vec *supp_pa_e = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec)); */
/**/
/*   for (unsigned int i = 0; i < pbwt->n_haps; i++) { */
/*     pa[i] = i; */
/*     da[i] = 0; */
/*     kv_init(supp_b[i]); */
/*     kv_init(supp_e[i]); */
/*     kv_init(supp_pa_b[i]); */
/*     kv_init(supp_pa_e[i]); */
/*     kv_init(supp_da_b[i]); */
/*   } */
/**/
/*   uint32_t c = 0; */
/**/
/*   uint8_t *c_col = (uint8_t *)malloc(pbwt->n_haps * sizeof(uint8_t)); */
/**/
/*   column_vec pbwt_cols; */
/*   kv_init(pbwt_cols); */
/**/
/*   kv_init(pbwt->cols); */
/*   while (bcf_read(fp, hdr, rec) >= 0) { */
/*     pbwt_col ds_col = {0}; */
/**/
/*     pbwt->n_sites++; */
/*     fprintf(stderr, "reading site = %d\r", c); */
/**/
/*     bcf_unpack(rec, BCF_UN_ALL); */
/*     uint32_t n_c = 0; */
/*     int32_t *gt_arr = NULL, ngt_arr = 0; */
/*     int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr); */
/*     ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr); */
/*     int max_ploidy = ngt / nsmpl; */
/**/
/*     for (i = 0; i < nsmpl; i++) { */
/*       int32_t *ptr = gt_arr + i * max_ploidy; */
/*       for (j = 0; j < max_ploidy; j++) { */
/*         if (ptr[j] == bcf_int32_vector_end) */
/*           break; */
/**/
/*         if (bcf_gt_is_missing(ptr[j])) */
/*           exit(-1); */
/**/
/*         uint32_t idx = i + i + j; */
/*         c_col[idx] = bcf_gt_allele(ptr[j]); */
/*       } */
/*     } */
/*     free(gt_arr); */
/*     pbwt_col col = build_col(c_col, pa, da, pbwt->n_haps); */
/*     for (size_t r = 0; r < col.p.n; r++) { */
/*       kv_push(uint32_t, supp_b[c_arr_get(&col.b_pa, r)], c); */
/*       kv_push(uint32_t, supp_e[c_arr_get(&col.e_pa, r)], c); */
/*       if (r == 0) { */
/*         kv_push(uint32_t, supp_pa_b[c_arr_get(&col.b_pa, r)], pbwt->n_haps);
 */
/*         kv_push(uint32_t, supp_da_b[c_arr_get(&col.b_pa, r)], 0); */
/*       } else { */
/*         kv_push(uint32_t, supp_pa_b[c_arr_get(&col.b_pa, r)], */
/*                 pa[c_arr_get(&col.p, r) - 1]); */
/*         kv_push(uint32_t, supp_da_b[c_arr_get(&col.b_pa, r)], */
/*                 da[c_arr_get(&col.p, r)]); */
/*       } */
/*       if (r == col.p.n - 1) { */
/*         kv_push(uint32_t, supp_pa_e[c_arr_get(&col.e_pa, r)], pbwt->n_haps);
 */
/*       } else { */
/*         kv_push(uint32_t, supp_pa_e[c_arr_get(&col.e_pa, r)], */
/*                 pa[c_arr_get(&col.p, r + 1)]); */
/*       } */
/*     } */
/**/
/*     memcpy(l_pa, pa, pbwt->n_haps * sizeof(uint32_t)); */
/*     memcpy(l_da, da, pbwt->n_haps * sizeof(uint32_t)); */
/**/
/*     KV_PUSH_COL_VEC(pbwt->cols, col); */
/**/
/*     pbwt_update(c_col, &pa, &da, pbwt->n_haps); */
/*     // free_pbwt_col(&col); */
/*     c++; */
/*   } */
/**/
/*   // pbwt->cols = pbwt_cols; */
/*   pbwt_col col = build_col(c_col, pa, da, pbwt->n_haps); */
/*   KV_PUSH_COL_VEC(pbwt->cols, col); */
/**/
/*   fprintf(stderr, "Building phi\n"); */
/*   memcpy(l_pa, pa, pbwt->n_haps * sizeof(uint32_t)); */
/*   memcpy(l_da, da, pbwt->n_haps * sizeof(uint32_t)); */
/*   for (size_t r = 0; r < col.p.n; r++) { */
/*     kv_push(uint32_t, supp_b[c_arr_get(&col.b_pa, r)], c); */
/*     kv_push(uint32_t, supp_e[c_arr_get(&col.e_pa, r)], c); */
/*     if (r == 0) { */
/*       kv_push(uint32_t, supp_pa_b[c_arr_get(&col.b_pa, r)], pbwt->n_haps); */
/*       kv_push(uint32_t, supp_da_b[c_arr_get(&col.b_pa, r)], 0); */
/*     } else { */
/*       kv_push(uint32_t, supp_pa_b[c_arr_get(&col.b_pa, r)], */
/*               pa[c_arr_get(&col.p, r) - 1]); */
/*       kv_push(uint32_t, supp_da_b[c_arr_get(&col.b_pa, r)], */
/*               da[c_arr_get(&col.p, r)]); */
/*     } */
/*     if (r == col.p.n - 1) { */
/*       kv_push(uint32_t, supp_pa_e[c_arr_get(&col.e_pa, r)], pbwt->n_haps); */
/*     } else { */
/*       kv_push(uint32_t, supp_pa_e[c_arr_get(&col.e_pa, r)], */
/*               pa[c_arr_get(&col.p, r + 1)]); */
/*     } */
/*   } */
/**/
/*   for (size_t i = 0; i < pbwt->n_haps; i++) { */
/*     if (i == 0) { */
/*       if (supp_pa_b[l_pa[i]].n == 0 || */
/*           supp_pa_b[l_pa[i]].a[supp_pa_b[l_pa[i]].n - 1] != pbwt->n_haps) {
 */
/*         kv_push(uint32_t, supp_pa_b[l_pa[i]], pbwt->n_haps); */
/*         kv_push(uint32_t, supp_da_b[l_pa[i]], 0); */
/*       } */
/*     } else { */
/*       if (supp_pa_b[l_pa[i]].n == 0 || */
/*           supp_pa_b[l_pa[i]].a[supp_pa_b[l_pa[i]].n - 1] != l_pa[i - 1]) { */
/*         kv_push(uint32_t, supp_pa_b[l_pa[i]], l_pa[i - 1]); */
/*         kv_push(uint32_t, supp_da_b[l_pa[i]], l_da[i]); */
/*       } */
/*     } */
/**/
/*     if (i == pbwt->n_haps - 1) { */
/*       if (supp_pa_e[l_pa[i]].n == 0 || */
/*           supp_pa_e[l_pa[i]].a[supp_pa_b[l_pa[i]].n - 1] != pbwt->n_haps) {
 */
/*         kv_push(uint32_t, supp_pa_e[l_pa[i]], pbwt->n_haps); */
/*       } */
/*     } else { */
/*       if (supp_pa_e[l_pa[i]].n == 0 || */
/*           supp_pa_e[l_pa[i]].a[supp_pa_e[l_pa[i]].n - 1] != l_pa[i + 1]) { */
/*         kv_push(uint32_t, supp_pa_e[l_pa[i]], l_pa[i + 1]); */
/*       } */
/*     } */
/*   } */
/**/
/*   build_phi(&pbwt->phi, supp_b, supp_e, supp_pa_b, supp_pa_e, supp_da_b, */
/*             pbwt->n_haps, pbwt->n_sites + 1); */
/**/
/*   for (size_t i = 0; i < pbwt->n_haps; i++) { */
/*     kv_destroy(supp_b[i]); */
/*     kv_destroy(supp_e[i]); */
/*     kv_destroy(supp_pa_b[i]); */
/*     kv_destroy(supp_pa_e[i]); */
/*     kv_destroy(supp_da_b[i]); */
/*   } */
/*   free(supp_b); */
/*   free(supp_e); */
/*   free(supp_da_b); */
/*   free(supp_pa_b); */
/*   free(supp_pa_e); */
/**/
/*   fprintf(stderr, "\n# sites = %d\n", pbwt->n_sites); */
/**/
/*   free(c_col); */
/*   free(pa); */
/*   free(da); */
/*   free(l_pa); */
/*   free(l_da); */
/**/
/*   bcf_destroy(rec); */
/*   bcf_hdr_destroy(hdr); */
/*   hts_close(fp); */
/* } */

void pbwt_build_af(pbwt *pbwt, std::string fref, ref_haplotype_set &H,
                   variant_map &V, std::string &region) {

  tac.clock();
  bool keep_mono = false;
  vrb.bullet("building mu-PBWT ");
  bcf_srs_t *sr = bcf_sr_init();
  sr->require_index = 1;
  if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1)
    vrb.error("Impossible to jump to region [" + region + "] in [" + fref +
              "]");
  if (bcf_sr_set_targets(sr, region.c_str(), 0, 0) == -1)
    vrb.error("Impossible to set target region [" + region + "] in [" + fref +
              "]");

  /* if (nthreads > 1) */
  /*   bcf_sr_set_threads(sr, nthreads); */

  if (!(bcf_sr_add_reader(sr, fref.c_str()))) {
    if (sr->errnum != idx_load_failed)
      vrb.error("Failed to open file: " + fref + "");
    else
      vrb.error("Failed to load index of the file: " + fref + "");
  }

  uint32_t n_ref_samples = bcf_hdr_nsamples(sr->readers[0].header);
  pbwt->n_haps = n_ref_samples * 2;
  pbwt->n_sites = 0;

  fprintf(stderr, "# haplotypes = %d\n", pbwt->n_haps);

  uint32_t *pa = (uint32_t *)malloc(pbwt->n_haps * sizeof(uint32_t));
  uint32_t *da = (uint32_t *)malloc(pbwt->n_haps * sizeof(uint32_t));
  uint32_t *l_pa = (uint32_t *)malloc(pbwt->n_haps * sizeof(uint32_t));
  uint32_t *l_da = (uint32_t *)malloc(pbwt->n_haps * sizeof(uint32_t));

  // uint32_t *thr = (uint32_t *)malloc(pbwt.n_haps * sizeof(uint32_t));
  // uint32_t *thr = NULL;

  int_vec *supp_b = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec));
  int_vec *supp_e = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec));

  int_vec *supp_pa_b = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec));
  int_vec *supp_da_b = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec));
  int_vec *supp_pa_e = (int_vec *)malloc(pbwt->n_haps * sizeof(int_vec));

  for (unsigned int i = 0; i < pbwt->n_haps; i++) {
    pa[i] = i;
    da[i] = 0;
    kv_init(supp_b[i]);
    kv_init(supp_e[i]);
    kv_init(supp_pa_b[i]);
    kv_init(supp_pa_e[i]);
    kv_init(supp_da_b[i]);
  }

  uint32_t c = 0;

  uint8_t *c_col = (uint8_t *)malloc(pbwt->n_haps * sizeof(uint8_t));

  column_vec pbwt_cols;
  kv_init(pbwt_cols);

  kv_init(pbwt->cols);

  unsigned int i_site = 0, i_common = 0, n_ref_unphased = 0;
  int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
  bcf1_t *line_ref;
  float prog_step = 1.0 / H.n_tot_sites;
  float prog_bar = 0.0;
  sr->max_unpack = BCF_UN_ALL;

  unsigned int cref = 0, calt = 0;
  int line_max_ploidy = 0;
  int idx_ref_hap = 0;
  bool a = false;
  int *ptr;
  float *ptr_f;

#ifdef __XSI__
  c_xcf *c_xcf_p = c_xcf_new();
  c_xcf_add_readers(c_xcf_p, sr);
#endif

  int rAC = 0, nAC = 0, *vAC = NULL;
  int rAN = 0, nAN = 0, *vAN = NULL;
  int ones = 0;
  int zeros = 0;
  auto skip1 = 0;
  auto skip2 = 0;

  std::vector<int> ploidy_ref_samples;

  while (bcf_sr_next_line(sr)) {
    if (bcf_sr_get_line(sr, 0)->n_allele != 2) {
      skip1++;
      continue; // always ref
    }
    line_ref = bcf_sr_get_line(sr, 0);

    rAC = bcf_get_info_int32(sr->readers[0].header, line_ref, "AC", &vAC, &nAC);
    rAN = bcf_get_info_int32(sr->readers[0].header, line_ref, "AN", &vAN, &nAN);
    if ((nAC != 1) || (nAN != 1))
      vrb.error(
          "VCF for reference panel needs AC/AN INFO fields to be present");
    calt = vAC[0];
    cref = (vAN[0] - vAC[0]);

    if (std::min(calt, cref) == 0 && !keep_mono) {
      skip2++;
      continue;
    }

#ifdef __XSI__
    ngt_ref = c_xcf_get_genotypes(c_xcf_p, 0, sr->readers[0].header, line_ref,
                                  &gt_arr_ref, &ngt_arr_ref);
#else
    ngt_ref = bcf_get_genotypes(sr->readers[0].header, line_ref, &gt_arr_ref,
                                &ngt_arr_ref);
#endif

    pbwt_col ds_col = {0};
    idx_ref_hap = 0;
    pbwt->n_sites++;

    if (i_site == 0) {
      int max_ploidy_ref = 0;
      int n_ref_haploid = 0;

      if (gt_arr_ref == nullptr || ngt_ref == 0)
        vrb.error("Error setting ploidy while reading the first record.");

      max_ploidy_ref = ngt_ref / n_ref_samples;
      if (max_ploidy_ref < 1 || max_ploidy_ref > 2 ||
          ngt_ref % n_ref_samples != 0)
        vrb.error("Max ploidy of the reference panel cannot be set to neither "
                  "1 or 2.");

      ploidy_ref_samples = std::vector<int>(n_ref_samples, 2);

      if (max_ploidy_ref == 1) {
        n_ref_haploid = n_ref_samples;
        for (int i = 0; i < n_ref_samples; ++i)
          ploidy_ref_samples[i] = 1;
      } else {
        for (int i = 0; i < n_ref_samples; ++i) {
          if (gt_arr_ref[max_ploidy_ref * i + 1] == bcf_int32_vector_end) {
            ploidy_ref_samples[i] = 1;
            ++n_ref_haploid;
          }
        }
      }
    }
    line_max_ploidy = ngt_ref / n_ref_samples;
    idx_ref_hap = 0, cref = 0, calt = 0;

    for (int i = 0; i < n_ref_samples; ++i) {
      ptr = gt_arr_ref + i * line_max_ploidy;
      for (int j = 0; j < ploidy_ref_samples[i]; j++) {

        ++idx_ref_hap;
        uint32_t idx = i + i + j;
        c_col[idx] = bcf_gt_allele(ptr[j]);
      }
    }

    pbwt_col col = build_col(c_col, pa, da, pbwt->n_haps);
    for (size_t r = 0; r < col.p.n; r++) {
      kv_push(uint32_t, supp_b[c_arr_get(&col.b_pa, r)], c);
      kv_push(uint32_t, supp_e[c_arr_get(&col.e_pa, r)], c);
      if (r == 0) {
        kv_push(uint32_t, supp_pa_b[c_arr_get(&col.b_pa, r)], pbwt->n_haps);
        kv_push(uint32_t, supp_da_b[c_arr_get(&col.b_pa, r)], 0);
      } else {
        kv_push(uint32_t, supp_pa_b[c_arr_get(&col.b_pa, r)],
                pa[c_arr_get(&col.p, r) - 1]);
        kv_push(uint32_t, supp_da_b[c_arr_get(&col.b_pa, r)],
                da[c_arr_get(&col.p, r)]);
      }
      if (r == col.p.n - 1) {
        kv_push(uint32_t, supp_pa_e[c_arr_get(&col.e_pa, r)], pbwt->n_haps);
      } else {
        kv_push(uint32_t, supp_pa_e[c_arr_get(&col.e_pa, r)],
                pa[c_arr_get(&col.p, r + 1)]);
      }
    }

    memcpy(l_pa, pa, pbwt->n_haps * sizeof(uint32_t));
    memcpy(l_da, da, pbwt->n_haps * sizeof(uint32_t));

    KV_PUSH_COL_VEC(pbwt->cols, col);

    pbwt_update(c_col, &pa, &da, pbwt->n_haps);
    // free_pbwt_col(&col);
    c++;
    /* if ((cref != V.vec_pos[i_site]->cref) || (calt !=
     * V.vec_pos[i_site]->calt)) */
    /*   vrb.error("AC/AN INFO fields in VCF are inconsistent with GT field, "
     */
    /*             "update the values in the VCF"); */
    i_site++;
    prog_bar += prog_step;
    vrb.progress(" * mu-PBWT building  ", prog_bar);
  }

  // pbwt->cols = pbwt_cols;
  pbwt_col col = build_col(c_col, pa, da, pbwt->n_haps);
  KV_PUSH_COL_VEC(pbwt->cols, col);

  fprintf(stderr, "Building phi\n");
  memcpy(l_pa, pa, pbwt->n_haps * sizeof(uint32_t));
  memcpy(l_da, da, pbwt->n_haps * sizeof(uint32_t));
  for (size_t r = 0; r < col.p.n; r++) {
    kv_push(uint32_t, supp_b[c_arr_get(&col.b_pa, r)], c);
    kv_push(uint32_t, supp_e[c_arr_get(&col.e_pa, r)], c);
    if (r == 0) {
      kv_push(uint32_t, supp_pa_b[c_arr_get(&col.b_pa, r)], pbwt->n_haps);
      kv_push(uint32_t, supp_da_b[c_arr_get(&col.b_pa, r)], 0);
    } else {
      kv_push(uint32_t, supp_pa_b[c_arr_get(&col.b_pa, r)],
              pa[c_arr_get(&col.p, r) - 1]);
      kv_push(uint32_t, supp_da_b[c_arr_get(&col.b_pa, r)],
              da[c_arr_get(&col.p, r)]);
    }
    if (r == col.p.n - 1) {
      kv_push(uint32_t, supp_pa_e[c_arr_get(&col.e_pa, r)], pbwt->n_haps);
    } else {
      kv_push(uint32_t, supp_pa_e[c_arr_get(&col.e_pa, r)],
              pa[c_arr_get(&col.p, r + 1)]);
    }
  }

  for (size_t i = 0; i < pbwt->n_haps; i++) {
    if (i == 0) {
      if (supp_pa_b[l_pa[i]].n == 0 ||
          supp_pa_b[l_pa[i]].a[supp_pa_b[l_pa[i]].n - 1] != pbwt->n_haps) {
        kv_push(uint32_t, supp_pa_b[l_pa[i]], pbwt->n_haps);
        kv_push(uint32_t, supp_da_b[l_pa[i]], 0);
      }
    } else {
      if (supp_pa_b[l_pa[i]].n == 0 ||
          supp_pa_b[l_pa[i]].a[supp_pa_b[l_pa[i]].n - 1] != l_pa[i - 1]) {
        kv_push(uint32_t, supp_pa_b[l_pa[i]], l_pa[i - 1]);
        kv_push(uint32_t, supp_da_b[l_pa[i]], l_da[i]);
      }
    }

    if (i == pbwt->n_haps - 1) {
      if (supp_pa_e[l_pa[i]].n == 0 ||
          supp_pa_e[l_pa[i]].a[supp_pa_b[l_pa[i]].n - 1] != pbwt->n_haps) {
        kv_push(uint32_t, supp_pa_e[l_pa[i]], pbwt->n_haps);
      }
    } else {
      if (supp_pa_e[l_pa[i]].n == 0 ||
          supp_pa_e[l_pa[i]].a[supp_pa_e[l_pa[i]].n - 1] != l_pa[i + 1]) {
        kv_push(uint32_t, supp_pa_e[l_pa[i]], l_pa[i + 1]);
      }
    }
  }

  build_phi(&pbwt->phid, supp_b, supp_e, supp_pa_b, supp_pa_e, supp_da_b,
            pbwt->n_haps, pbwt->n_sites + 1);

  for (size_t i = 0; i < pbwt->n_haps; i++) {
    kv_destroy(supp_b[i]);
    kv_destroy(supp_e[i]);
    kv_destroy(supp_pa_b[i]);
    kv_destroy(supp_pa_e[i]);
    kv_destroy(supp_da_b[i]);
  }
  free(supp_b);
  free(supp_e);
  free(supp_da_b);
  free(supp_pa_b);
  free(supp_pa_e);

  fprintf(stderr, "\n# sites = %d\n", pbwt->n_sites);

  free(c_col);
  free(pa);
  free(da);
  free(l_pa);
  free(l_da);

  /* bcf_destroy(rec); */
  /* bcf_hdr_destroy(hdr); */
  /* hts_close(fp); */

#ifdef __XSI__
  c_xcf_delete(c_xcf_p);
#endif
  free(gt_arr_ref);
  vrb.bullet("mu-PBWT building (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) +
             "s)");
}

void pbwt_update(uint8_t *col, uint32_t **pa, uint32_t **da, uint32_t n_h) {
  uint32_t *n_pa = (uint32_t *)malloc(n_h * sizeof(uint32_t));
  uint32_t *n_da = (uint32_t *)malloc(n_h * sizeof(uint32_t));

  uint32_t c_z = 0;
  uint32_t lcs = -1;

  for (size_t i = 0; i < n_h; i++) {
    lcs = U32_MIN(lcs, (*da)[i]);
    if (col[(*pa)[i]] == 0) {
      n_pa[c_z] = (*pa)[i];
      n_da[c_z] = lcs + 1;
      c_z++;
      lcs = UINT32_MAX;
    }
  }
  uint32_t c_o = 0;
  lcs = -1;
  for (size_t i = 0; i < n_h; i++) {
    lcs = U32_MIN(lcs, (*da)[i]);
    if (col[(*pa)[i]] == 1) {
      n_pa[c_z + c_o] = (*pa)[i];
      n_da[c_z + c_o] = lcs + 1;
      c_o++;
      lcs = UINT32_MAX;
    }
  }
  n_da[0] = 0;
  if (c_z != n_h) {
    n_da[c_z] = 0;
  }
  /* DBG_PRINT_ARRAY("%d", n_pa, n_h); */
  /* DBG_PRINT_ARRAY("%d", n_da, n_h); */
  free(*pa);
  free(*da);
  *da = n_da;
  *pa = n_pa;
}

void get_pbwt_col(uint8_t *col, uint32_t **pa, uint8_t **p_col, uint32_t n_h) {
  for (size_t i = 0; i < n_h; i++) {
    (*p_col)[i] = col[(*pa)[i]];
  }
}

uv_res get_uv(const pbwt_col *col, uint32_t r) {
  if (r == 0)
    return (uv_res){0, 0};

  uint32_t a = c_arr_get(&col->uv, r);
  uint32_t b = c_arr_get(&col->uv, r - 1);

  uint32_t u, v;
  if (r & 1) {
    u = a;
    v = b;
  } else {
    u = b;
    v = a;
  }

  if (!col->zero) {
    return (uv_res){v, u};
  }

  return (uv_res){u, v};
}

// assuming s is always the current run symboil, hence s is removed
uint32_t fl(const pbwt_col *col, uint32_t i, uint32_t r) {
  uv_res uv = get_uv(col, r);
  uint32_t run_start = c_arr_get(&col->p, r);
  uint32_t offset = i - run_start;
  uint8_t s = get_ns(col->zero, r);
  if (s == 0) {
    return uv.u + offset;
  } else {
    return col->n_zeros + uv.v + offset;
  }
}

uint32_t lf(const column_vec *cols, uint32_t col_idx, uint32_t i) {
  if (col_idx == 0)
    return 0;

  const pbwt_col *prev_col = &cols->a[col_idx - 1];
  uint32_t c0 = prev_col->n_zeros;

  int is_one = (i >= c0);
  uint32_t target_rank = is_one ? (i - c0) : i;

  size_t pos = 0;
  size_t n = prev_col->p.n;
  if (n > 0) {
    size_t step = 1ULL << (31 - __builtin_clz((uint32_t)n));
    for (; step > 0; step >>= 1) {
      size_t next_pos = pos + step;
      if (next_pos < n) {
        uv_res uv = get_uv(prev_col, (uint32_t)next_pos);
        uint32_t val = is_one ? uv.v : uv.u;
        if (val <= target_rank) {
          pos = next_pos;
        }
      }
    }
  }

  uv_res final_uv = get_uv(prev_col, (uint32_t)pos);
  uint32_t base_rank = is_one ? final_uv.v : final_uv.u;
  uint32_t offset = target_rank - base_rank;

  return c_arr_get(&prev_col->p, (uint32_t)pos) + offset;
}

uint32_t get_l(const pbwt *pbwt, const uint8_t *q, uint32_t i, uint32_t c_i) {
  int32_t t_i = (int32_t)i;
  uint32_t t_lf = c_i;
  uint32_t t_l = 0;
  uint32_t last_r = 0;

  const pbwt_col *all_cols = pbwt->cols.a;

  while (t_i >= 0) {
    __builtin_prefetch(&all_cols[t_i - 1], 0, 3);
    const pbwt_col *col = &all_cols[t_i];

    uint32_t r = get_r_with_hint(col, t_lf, last_r);
    last_r = r;

    if (q[t_i] != get_ns(col->zero, r)) {
      break;
    }

    if (t_i > 0) {
      t_lf = lf(&pbwt->cols, (uint32_t)t_i, t_lf);
    }

    t_i--;
    t_l++;
  }

  return t_l;
}
int_vec get_haps_ms(const pbwt *pbwt, uint32_t p, uint32_t l, uint32_t c) {
  int_vec res;
  kv_init(res);

  kv_push(uint32_t, res, p);
  uint32_t t_p = p;

  bool d = true;
  uint32_t d_p = 0;
  bool u = true;
  uint32_t u_p = 0;

  /* map_res d_m = pbwt_for(pbwt, c_i, c_r, c); */
  /* d_m.i += 1; */
  /* map_res u_m = pbwt_for(pbwt, c_i, c_r, c); */

  while (d) {
    d_p = phi_inv_f(&pbwt->phid, t_p, c + 1);
    DBG("testing down %d", d_p);
    if (d_p == pbwt->n_haps)
      break;
    if (phi_l(&pbwt->phid, d_p, c + 1) >= l) {
      kv_push(uint32_t, res, d_p);
      t_p = d_p;
    } else {
      d = false;
    }
  }

  t_p = p;

  while (u) {
    u_p = phi_f(&pbwt->phid, t_p, c + 1);
    DBG("testing up %d", u_p);

    if (u_p == pbwt->n_haps)
      break;
    if (phi_l(&pbwt->phid, t_p, c + 1) >= l) {
      kv_push(uint32_t, res, u_p);
      t_p = u_p;
    } else {
      u = false;
    }
  }

  // sort_int_vec(&res);
  return res;
}
uint32_t get_l_u(const pbwt *pbwt, uint32_t p, uint32_t n, uint32_t c) {
  if (p >= pbwt->n_haps || n == 0)
    return 0;

  uint32_t m = UINT32_MAX;
  uint32_t t_p = p;
  const uint32_t sentinel = pbwt->n_haps;

  for (uint32_t i = 0; i < n; i++) {
    uint32_t t_l = phi_l(&pbwt->phid, t_p, c);
    if (t_l < m)
      m = t_l;
    if (m == 0)
      return 0;

    t_p = phi_f(&pbwt->phid, t_p, c);

    if (t_p >= sentinel)
      break;
  }
  return m;
}
uint32_t get_l_d(const pbwt *pbwt, uint32_t p, uint32_t n, uint32_t c) {
  uint32_t m = UINT32_MAX;
  uint32_t t_p = p;
  uint32_t sentinel = pbwt->n_haps;

  for (size_t i = 0; i < n; i++) {
    uint32_t next_p = phi_inv_f(&pbwt->phid, t_p, c);

    if (next_p >= sentinel)
      break;

    uint32_t t_l = phi_l(&pbwt->phid, next_p, c);

    if (t_l < m)
      m = t_l;
    if (m == 0)
      break;

    t_p = next_p;
  }
  return (m == UINT32_MAX) ? 0 : m;
}
match_vec compute_smem(const pbwt *pbwt, const uint8_t *q, uint32_t q_n) {
  match_vec res;
  kv_init(res);

  uint32_t p = 0, p_p = 0;
  uint32_t l = 0, p_l = 0;

  pbwt_col *all_cols = pbwt->cols.a;
  uint32_t n_haps = pbwt->n_haps;
  uint32_t n_sites = pbwt->n_sites;

  pbwt_col *c_col = &all_cols[0];
  uint32_t c_p = c_arr_get(&c_col->e_pa, c_col->e_pa.n - 1);
  uint32_t c_i = n_haps - 1;
  uint32_t c_r = c_col->e_pa.n - 1;
  uint8_t c_s = get_ns(c_col->zero, c_r);
  uint32_t p_i = c_i;

  for (size_t i = 0; i < n_sites; i++) {
    c_col = &all_cols[i];
    pbwt_col *next_col = (i < n_sites - 1) ? &all_cols[i + 1] : NULL;

    uint8_t q_s = q[i];

    if (LIKELY(q_s == c_s)) {
      p = c_p;
      l = (i != 0 && p_p == p) ? p_l + 1 : 1;
    } else {
      uint32_t t = c_arr_get(&c_col->t, c_r);
      uint32_t p_n = c_col->p.n;

      if (UNLIKELY(p_n == 1)) {
        p = n_haps;
        l = 0;
        if (next_col) {
          c_p = c_arr_get(&next_col->e_pa, next_col->e_pa.n - 1);
          c_i = n_haps - 1;
          c_r = next_col->e_pa.n - 1;
          c_s = get_ns(next_col->zero, c_r);
        }
        goto check_match;
      }

      uint32_t d;
      bool is_up = (c_r != 0 && ((c_i < t) || (c_r == p_n - 1)));

      if (is_up) {
        p_i = c_i;
        c_i = c_arr_get(&c_col->p, c_r) - 1;
        c_p = c_arr_get(&c_col->e_pa, c_r - 1);
        c_r--;
        d = p_i - c_i;
      } else {
        p_i = c_i;
        c_i = c_arr_get(&c_col->p, c_r + 1);
        c_p = c_arr_get(&c_col->b_pa, c_r + 1);
        c_r++;
        d = c_i - p_i;
      }

      if (d < 90) {
        uint32_t m =
            is_up ? get_l_u(pbwt, p_p, d, i) : get_l_d(pbwt, p_p, d, i);
        l = MIN(m, p_l) + 1;
      } else {
        l = get_l(pbwt, q, i, c_i);
      }

      p = c_p;
    }

    if (next_col) {
      c_i = fl(c_col, c_i, c_r);
      c_r = get_r(next_col, c_i);
      c_s = get_ns(next_col->zero, c_r);
    }

  check_match:
    if (i > 0 && p_p != p && p_l > 0 && p_l >= l) {
      match m = {(uint32_t)(i - 1), (uint32_t)p_l, q_n,
                 get_haps_ms(pbwt, p_p, p_l, i - 1)};
      KV_PUSH_MATCH_VEC(res, m);
    }

    if (i == n_sites - 1 && l != 0) {
      match m = {(uint32_t)i, (uint32_t)l, q_n, get_haps_ms(pbwt, p, l, i)};
      KV_PUSH_MATCH_VEC(res, m);
    }

    p_p = p;
    p_l = l;
  }
  return res;
}

void pbwt_query(const char *pbwt_filename, const char *filename, int threads) {
  omp_set_num_threads(threads);
  clock_t START = clock();
  struct timespec STARTM, ENDM;
  clock_gettime(CLOCK_MONOTONIC, &STARTM);

  pbwt pbwt;
  memset(&pbwt, 0, sizeof(pbwt));

  FILE *fpb = fopen(pbwt_filename, "rb");
  if (!fpb) {
    exit(-1);
  }
  setvbuf(fpb, NULL, _IOFBF, 1024 * 1024 * 4);

  pbwt_deserialize(fpb, &pbwt);
  fclose(fpb);
  clock_gettime(CLOCK_MONOTONIC, &ENDM);

  double wall_time =
      ENDM.tv_sec - STARTM.tv_sec + (ENDM.tv_nsec - STARTM.tv_nsec) / 1e9;
  fprintf(stderr, "MuPBWT loaded in %fs\n", wall_time);

  htsFile *fp = hts_open(filename, "rb");

  if (fp == NULL) {
    fprintf(stderr, "%s NOT FOUND\n", filename);
    exit(-1);
  }
  clock_gettime(CLOCK_MONOTONIC, &STARTM);

  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init();

  uint32_t n = bcf_hdr_nsamples(hdr) * 2;
  fprintf(stderr, "Found %u queries\n", n);

  uint32_t check_s = 0;
  uint32_t n_sites = pbwt.n_sites;

  uint8_t *q_p = NULL;
  size_t mat_size = (size_t)n * n_sites * sizeof(uint8_t);

  if (posix_memalign((void **)&q_p, 64, mat_size) != 0) {
    fprintf(stderr, "Memory allocation error for q_p\n");
    exit(-1);
  }
  memset(q_p, 0, mat_size);

  int32_t *gt_arr = NULL;
  int32_t ngt_arr = 0;

  while (bcf_read(fp, hdr, rec) >= 0 && check_s < n_sites) {
    bcf_unpack(rec, BCF_UN_FMT);

    int nsmpl = bcf_hdr_nsamples(hdr);
    int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
    int max_ploidy = ngt / nsmpl;

    for (int idx = 0; idx < nsmpl; idx++) {
      int32_t *ptr = gt_arr + idx * max_ploidy;
      for (int j = 0; j < max_ploidy; j++) {
        if (ptr[j] == bcf_int32_vector_end)
          break;
        if (bcf_gt_is_missing(ptr[j]))
          exit(-1);

        if (bcf_gt_allele(ptr[j]) == 1) {
          uint32_t q_idx = idx * 2 + j;
          q_p[(size_t)q_idx * n_sites + check_s] = 1;
        }
      }
    }
    check_s++;
  }

  free(gt_arr);
  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  bcf_close(fp);

  if (check_s != n_sites) {
    fprintf(stderr, "ERROR: pbwt sites = %u and query sites = %u\n", n_sites,
            check_s);
    exit(-1);
  }
  clock_gettime(CLOCK_MONOTONIC, &ENDM);

  wall_time =
      ENDM.tv_sec - STARTM.tv_sec + (ENDM.tv_nsec - STARTM.tv_nsec) / 1e9;
  fprintf(stderr, "Queries loaded in %fs\n", wall_time);
  clock_gettime(CLOCK_MONOTONIC, &STARTM);

#pragma omp parallel for default(none) shared(pbwt, q_p, n, n_sites, stdout)   \
    schedule(dynamic) ordered
  for (size_t i = 0; i < n; i++) {
    match_vec local_match = compute_smem(&pbwt, &q_p[(size_t)i * n_sites], i);

#pragma omp ordered
    {
      match_vec_print(&local_match);
    }

    match_vec_free(&local_match);
  }
  clock_gettime(CLOCK_MONOTONIC, &ENDM);

  wall_time =
      ENDM.tv_sec - STARTM.tv_sec + (ENDM.tv_nsec - STARTM.tv_nsec) / 1e9;
  fprintf(stderr, "Smems computed in %fs\n", wall_time);
  free_pbwt(&pbwt);
  free(q_p);
}

void free_pbwt(pbwt *p) {
  if (!p)
    return;

  for (size_t i = 0; i < p->n_sites + 1; i++) {
    free_pbwt_col(&p->cols.a[i]);
  }
  kv_destroy(p->cols);

  phi_free(&p->phid);
  p->n_haps = 0;
  p->n_sites = 0;
}

int pbwt_serialize_gz(const char *path, const pbwt *p) {
  if (!path || !p)
    return -1;

  gzFile fp = gzopen(path, "wb");
  if (fp == NULL)
    return -1;

  gzwrite(fp, &p->n_haps, sizeof(uint32_t));
  gzwrite(fp, &p->n_sites, sizeof(uint32_t));

  phi_serialize_gz(&p->phid, fp);
  for (size_t i = 0; i < p->cols.n; i++) {
    pbwt_col_serialize_gz(fp, &p->cols.a[i]);
  }

  gzclose(fp);
  return 0;
}

int pbwt_deserialize_gz(const char *path, pbwt *p) {
  if (!path || !p)
    return -1;

  gzFile fp = gzopen(path, "rb");
  if (fp == NULL)
    return -1;

  if (gzread(fp, &p->n_haps, sizeof(uint32_t)) != sizeof(uint32_t))
    goto error;
  if (gzread(fp, &p->n_sites, sizeof(uint32_t)) != sizeof(uint32_t))
    goto error;
  if (phi_deserialize_gz(&p->phid, fp) != 0)
    goto error;

  p->cols.n = p->n_sites + 1;
  p->cols.m = p->n_sites + 1;
  p->cols.a = (pbwt_col *)malloc(sizeof(pbwt_col) * p->cols.n);

  for (size_t i = 0; i < p->cols.n; i++) {
    if (pbwt_col_deserialize_gz(fp, &p->cols.a[i]) != 0)
      goto error;
  }

  gzclose(fp);
  return 0;

error:
  gzclose(fp);
  return -1;
}

int pbwt_serialize(FILE *fp, const pbwt *p) {
  if (!fp || !p)
    return -1;
  fwrite(&p->n_haps, sizeof(uint32_t), 1, fp);
  fwrite(&p->n_sites, sizeof(uint32_t), 1, fp);
  phi_serialize(&p->phid, fp);

  for (size_t i = 0; i < p->cols.n; i++) {
    pbwt_col_serialize(fp, &p->cols.a[i]);
  }
  return 0;
}

int pbwt_deserialize(FILE *fp, pbwt *p) {
  if (!fp || !p)
    return -1;
  fread(&p->n_haps, sizeof(uint32_t), 1, fp);
  fread(&p->n_sites, sizeof(uint32_t), 1, fp);
  p->cols.n = p->n_sites + 1;
  p->cols.m = p->n_sites + 1;
  p->cols.a = (pbwt_col *)malloc(sizeof(pbwt_col) * (p->n_sites + 1));
  phi_deserialize(&p->phid, fp);

  for (size_t i = 0; i < p->cols.n; i++) {
    // printf("column %zu\n", i);
    pbwt_col_deserialize(fp, &p->cols.a[i]);
  }
  return 0;
}

void print_pbwt(const pbwt *p) {
  if (!p)
    return;

  printf("Number of haplotypes: %u\n", p->n_haps);
  printf("Number of sites: %u\n", p->n_sites);

  for (size_t i = 0; i < p->cols.n; ++i) {
    printf("\n--- Column %zu ---\n", i);
    print_pbwt_col(&p->cols.a[i], p->n_haps);
  }
  phi_print(&p->phid);
}

void print_pa(const pbwt *p, uint32_t c) {
  int_vec pa;
  kv_init(pa);

  uint32_t s = c_arr_get(&p->cols.a[c].b_pa, 0);
  DBG("considering %d", s);
  kv_push(uint32_t, pa, s);
  uint32_t n = phi_inv_f(&p->phid, s, c);
  DBG("considering %d", n);
  while (n != p->n_haps) {
    kv_push(uint32_t, pa, n);
    n = phi_inv_f(&p->phid, n, c);
    DBG("considering %d", n);
  }
  print_kvec_uint32(&pa);
}

void print_da(const pbwt *p, uint32_t c) {
  int_vec da;
  kv_init(da);

  uint32_t s = c_arr_get(&p->cols.a[c].b_pa, 0);
  kv_push(uint32_t, da, phi_l(&p->phid, s, c));
  uint32_t n = phi_inv_f(&p->phid, s, c);
  while (n != p->n_haps) {
    kv_push(uint32_t, da, phi_l(&p->phid, n, c));
    n = phi_inv_f(&p->phid, n, c);
  }
  print_kvec_uint32(&da);
}

/* static int compare_groups(const void *a, const void *b) { */
/*   const smem_group *ga = (const smem_group *)a; */
/*   const smem_group *gb = (const smem_group *)b; */
/**/
/*   if (ga->p_p != gb->p_p) */
/*     return ga->p_p < gb->p_p ? -1 : 1; */
/*   if (ga->p_l != gb->p_l) */
/*     return ga->p_l < gb->p_l ? -1 : 1; */
/*   if (ga->c_i != gb->c_i) */
/*     return ga->c_i < gb->c_i ? -1 : 1; */
/*   if (ga->c_r != gb->c_r) */
/*     return ga->c_r < gb->c_r ? -1 : 1; */
/*   if (ga->c_p != gb->c_p) */
/*     return ga->c_p < gb->c_p ? -1 : 1; */
/*   if (ga->c_s != gb->c_s) */
/*     return ga->c_s < gb->c_s ? -1 : 1; */
/*   return 0; */
/* } */
/**/
/* static void emit_smems(match_vec *matches, const uint64_t *mask, */
/*                        uint32_t mask_words, uint32_t p_p, uint32_t p_l, */
/*                        uint32_t e, const pbwt *pbwt) { */
/*   if (p_p >= pbwt->n_haps) */
/*     return; */
/**/
/*   int_vec haps = get_haps_ms(pbwt, p_p, p_l, e); */
/**/
/*   for (uint32_t w = 0; w < mask_words; w++) { */
/*     uint64_t m = mask[w]; */
/*     while (m) { */
/*       int b = __builtin_ctzll(m); */
/*       uint32_t q_idx = w * 64 + b; */
/**/
/*       match match_obj = {.l = p_l, .q_n = q_idx, .e = e}; */
/*       kv_init(match_obj.h); */
/*       for (size_t hi = 0; hi < kv_size(haps); hi++) { */
/*         kv_push(uint32_t, match_obj.h, kv_A(haps, hi)); */
/*       } */
/*       KV_PUSH_MATCH_VEC(matches[q_idx], match_obj); */
/**/
/*       m &= m - 1; */
/*     } */
/*   } */
/*   kv_destroy(haps); */
/* } */
/* static void merge_or_push(group_vec *vec, smem_group gs, uint32_t mask_words)
 * { */
/*   for (size_t i = 0; i < kv_size(*vec); i++) { */
/*     smem_group *existing = &kv_A(*vec, i); */
/*     if (existing->p_p == gs.p_p && existing->p_l == gs.p_l && */
/*         existing->c_i == gs.c_i && existing->c_r == gs.c_r && */
/*         existing->c_p == gs.c_p && existing->c_s == gs.c_s) { */
/**/
/*       for (uint32_t w = 0; w < mask_words; w++) { */
/*         existing->q_mask[w] |= gs.q_mask[w]; */
/*       } */
/*       free(gs.q_mask); */
/*       return; */
/*     } */
/*   } */
/*   kv_push(smem_group, *vec, gs); */
/* } */
/**/
/* static inline void read_bcf_row_safe(htsFile *fp, bcf_hdr_t *hdr, bcf1_t
 * *rec, */
/*                                      int32_t **gt_arr, int32_t *ngt_arr, */
/*                                      uint8_t *q_row, uint32_t n) { */
/*   if (bcf_read(fp, hdr, rec) < 0) { */
/*     fprintf(stderr, "ERROR: VCF ended prematurely!\n"); */
/*     exit(-1); */
/*   } */
/*   bcf_unpack(rec, BCF_UN_FMT); */
/*   memset(q_row, 0, n); */
/**/
/*   int nsmpl = bcf_hdr_nsamples(hdr); */
/*   int ngt = bcf_get_genotypes(hdr, rec, gt_arr, ngt_arr); */
/*   int max_ploidy = ngt / nsmpl; */
/**/
/*   for (int idx = 0; idx < nsmpl; idx++) { */
/*     int32_t *ptr = (*gt_arr) + idx * max_ploidy; */
/*     for (int j = 0; j < max_ploidy; j++) { */
/*       if (ptr[j] == bcf_int32_vector_end) */
/*         break; */
/*       if (bcf_gt_is_missing(ptr[j])) */
/*         exit(-1); */
/*       if (bcf_gt_allele(ptr[j]) == 1) { */
/*         q_row[idx * 2 + j] = 1; */
/*       } */
/*     } */
/*   } */
/* } */

/* void pbwt_query_mq(const char *pbwt_filename, const char *filename, */
/*                    int threads) { */
/*   omp_set_num_threads(threads); */
/*   pbwt p; */
/*   memset(&p, 0, sizeof(pbwt)); */
/**/
/*   FILE *fpb = fopen(pbwt_filename, "rb"); */
/*   if (!fpb) { */
/*     fprintf(stderr, "Error loading PBWT file\n"); */
/*     exit(-1); */
/*   } */
/*   setvbuf(fpb, NULL, _IOFBF, 1024 * 1024 * 4); */
/**/
/*   fread(&p.n_haps, sizeof(uint32_t), 1, fpb); */
/*   fread(&p.n_sites, sizeof(uint32_t), 1, fpb); */
/*   phi_deserialize(&p.phid, fpb); */
/**/
/*   htsFile *fp = hts_open(filename, "rb"); */
/*   if (fp == NULL) { */
/*     fprintf(stderr, "%s NOT FOUND\n", filename); */
/*     exit(-1); */
/*   } */
/*   bcf_hdr_t *hdr = bcf_hdr_read(fp); */
/*   bcf1_t *rec = bcf_init(); */
/**/
/*   uint32_t n = bcf_hdr_nsamples(hdr) * 2; */
/*   fprintf(stderr, "Found %u queries\n", n); */
/**/
/*   match_vec *matches = (match_vec *)malloc(n * sizeof(match_vec)); */
/*   for (size_t q = 0; q < n; q++) { */
/*     kv_init(matches[q]); */
/*     kv_resize(match, matches[q], 256); */
/*   } */
/**/
/*   uint8_t *q_buffer[3]; */
/*   for (int k = 0; k < 3; k++) { */
/*     q_buffer[k] = (uint8_t *)malloc(n * sizeof(uint8_t)); */
/*   } */
/**/
/*   int32_t *gt_arr = NULL; */
/*   int32_t ngt_arr = 0; */
/**/
/*   if (p.n_sites > 0) */
/*     read_bcf_row_safe(fp, hdr, rec, &gt_arr, &ngt_arr, q_buffer[0], n); */
/*   if (p.n_sites > 1) */
/*     read_bcf_row_safe(fp, hdr, rec, &gt_arr, &ngt_arr, q_buffer[1], n); */
/**/
/*   pbwt_col cols[3]; */
/*   memset(cols, 0, sizeof(pbwt_col) * 3); */
/*   pbwt_col *c_col = &cols[0]; */
/*   pbwt_col *next_col = &cols[1]; */
/*   pbwt_col *next_next_col = &cols[2]; */
/**/
/*   pbwt_col_deserialize(fpb, c_col); */
/*   if (p.n_sites > 1) */
/*     pbwt_col_deserialize(fpb, next_col); */
/**/
/*   q_state *states = (q_state *)malloc(n * sizeof(q_state)); */
/*   uint32_t init_c_i = p.n_haps - 1; */
/*   uint32_t init_c_r = c_col->e_pa.n - 1; */
/*   uint32_t init_c_p = c_arr_get(&c_col->e_pa, init_c_r); */
/*   uint8_t init_c_s = get_ns(c_col->zero, init_c_r); */
/**/
/*   for (size_t q = 0; q < n; q++) { */
/*     states[q].p_p = 0; */
/*     states[q].p_l = 0; */
/*     states[q].c_i = init_c_i; */
/*     states[q].c_r = init_c_r; */
/*     states[q].c_p = init_c_p; */
/*     states[q].c_s = init_c_s; */
/*   } */
/**/
/*   uint32_t i = 0; */
/**/
/* #pragma omp parallel default(none) \ */
/*     shared(n, states, q_buffer, c_col, next_col, next_next_col, p, i,
 * matches, \ */
/*                fpb, fp, hdr, rec, gt_arr, ngt_arr) */
/*   { */
/*     while (i < p.n_sites) { */
/*       bool has_next = (i < p.n_sites - 1); */
/*       bool has_next_next = (i < p.n_sites - 2); */
/**/
/*       uint8_t *q_row = q_buffer[i % 3]; */
/**/
/* #pragma omp single nowait */
/*       { */
/*         if (has_next_next) { */
/*           pbwt_col_deserialize(fpb, next_next_col); */
/*           read_bcf_row_safe(fp, hdr, rec, &gt_arr, &ngt_arr, */
/*                             q_buffer[(i + 2) % 3], n); */
/*         } */
/*       } */
/**/
/* #pragma omp for schedule(dynamic, 64) */
/*       for (size_t q = 0; q < n; q++) { */
/*         uint8_t q_s = q_row[q]; */
/*         q_state s = states[q]; */
/*         uint32_t cur_p, cur_l; */
/**/
/*         if (LIKELY(q_s == s.c_s)) { */
/*           cur_p = s.c_p; */
/*           cur_l = (i != 0 && s.p_p == cur_p) ? s.p_l + 1 : 1; */
/*         } else { */
/*           uint32_t t = c_arr_get(&c_col->t, s.c_r); */
/*           uint32_t p_n = c_col->p.n; */
/**/
/*           if (UNLIKELY(p_n == 1)) { */
/*             cur_p = p.n_haps; */
/*             cur_l = 0; */
/*             if (has_next) { */
/*               s.c_p = c_arr_get(&next_col->e_pa, next_col->e_pa.n - 1); */
/*               s.c_i = p.n_haps - 1; */
/*               s.c_r = next_col->e_pa.n - 1; */
/*               s.c_s = get_ns(next_col->zero, s.c_r); */
/*             } */
/*             goto check_match; */
/*           } */
/**/
/*           uint32_t d; */
/*           bool is_up = (s.c_r != 0 && ((s.c_i < t) || (s.c_r == p_n - 1)));
 */
/*           uint32_t p_i = s.c_i; */
/**/
/*           if (is_up) { */
/*             s.c_i = c_arr_get(&c_col->p, s.c_r) - 1; */
/*             s.c_p = c_arr_get(&c_col->e_pa, s.c_r - 1); */
/*             s.c_r--; */
/*             d = p_i - s.c_i; */
/*           } else { */
/*             s.c_i = c_arr_get(&c_col->p, s.c_r + 1); */
/*             s.c_p = c_arr_get(&c_col->b_pa, s.c_r + 1); */
/*             s.c_r++; */
/*             d = s.c_i - p_i; */
/*           } */
/**/
/*           if (s.p_p < p.n_haps) { */
/*             uint32_t m = */
/*                 is_up ? get_l_u(&p, s.p_p, d, i) : get_l_d(&p, s.p_p, d, i);
 */
/*             cur_l = MIN(m, s.p_l) + 1; */
/*           } else { */
/*             cur_l = 1; */
/*           } */
/*           cur_p = s.c_p; */
/*         } */
/**/
/*         if (has_next) { */
/*           s.c_i = fl(c_col, s.c_i, s.c_r); */
/*           s.c_r = get_r(next_col, s.c_i); */
/*           s.c_s = get_ns(next_col->zero, s.c_r); */
/*         } */
/**/
/*       check_match: */
/*         if (i > 0 && s.p_p != cur_p && s.p_l > 0 && s.p_l >= cur_l) { */
/*           match m = {.l = s.p_l, */
/*                      .q_n = (uint32_t)q, */
/*                      .e = i - 1, */
/*                      .h = get_haps_ms(&p, s.p_p, s.p_l, i - 1)}; */
/*           KV_PUSH_MATCH_VEC(matches[q], m); */
/*         } */
/**/
/*         if (i == p.n_sites - 1 && cur_l != 0) { */
/*           match m = {.l = cur_l, */
/*                      .q_n = (uint32_t)q, */
/*                      .e = i, */
/*                      .h = get_haps_ms(&p, cur_p, cur_l, i)}; */
/*           KV_PUSH_MATCH_VEC(matches[q], m); */
/*         } */
/**/
/*         s.p_p = cur_p; */
/*         s.p_l = cur_l; */
/*         states[q] = s; */
/*       } */
/**/
/* #pragma omp barrier */
/**/
/* #pragma omp single */
/*       { */
/*         free_pbwt_col(c_col); */
/**/
/*         pbwt_col *temp = c_col; */
/*         c_col = next_col; */
/*         next_col = next_next_col; */
/*         next_next_col = temp; */
/**/
/*         i++; */
/*       } */
/*     } */
/*   } */
/**/
/*   if (c_col) */
/*     free_pbwt_col(c_col); */
/*   if (next_col && i <= p.n_sites) */
/*     free_pbwt_col(next_col); */
/**/
/*   free(states); */
/*   for (int k = 0; k < 3; k++) */
/*     free(q_buffer[k]); */
/*   free(gt_arr); */
/**/
/*   for (size_t q = 0; q < n; q++) { */
/*     match_vec_print(&matches[q]); */
/*     match_vec_free(&matches[q]); */
/*   } */
/*   free(matches); */
/**/
/*   bcf_destroy(rec); */
/*   bcf_hdr_destroy(hdr); */
/*   bcf_close(fp); */
/*   phi_free(&p.phid); */
/*   fclose(fpb); */
/* } */
/**/
/* void pbwt_query_m(const char *pbwt_filename, const char *filename, */
/*                   int threads) { */
/*   omp_set_num_threads(threads); */
/*   pbwt p; */
/*   memset(&p, 0, sizeof(pbwt)); */
/**/
/*   FILE *fpb = fopen(pbwt_filename, "rb"); */
/*   if (!fpb) { */
/*     fprintf(stderr, "Error loading PBWT file\n"); */
/*     exit(-1); */
/*   } */
/*   setvbuf(fpb, NULL, _IOFBF, 1024 * 1024 * 4); */
/**/
/*   fread(&p.n_haps, sizeof(uint32_t), 1, fpb); */
/*   fread(&p.n_sites, sizeof(uint32_t), 1, fpb); */
/*   phi_deserialize(&p.phid, fpb); */
/**/
/*   htsFile *fp = hts_open(filename, "rb"); */
/*   if (fp == NULL) { */
/*     fprintf(stderr, "%s NOT FOUND\n", filename); */
/*     exit(-1); */
/*   } */
/*   bcf_hdr_t *hdr = bcf_hdr_read(fp); */
/*   bcf1_t *rec = bcf_init(); */
/**/
/*   uint32_t n = bcf_hdr_nsamples(hdr) * 2; */
/*   fprintf(stderr, "Found %u queries\n", n); */
/**/
/*   match_vec *matches = (match_vec *)malloc(n * sizeof(match_vec)); */
/*   for (size_t q = 0; q < n; q++) { */
/*     kv_init(matches[q]); */
/*     kv_resize(match, matches[q], 256); */
/*   } */
/**/
/*   uint8_t *q_matrix = NULL; */
/*   size_t mat_size = (size_t)n * p.n_sites * sizeof(uint8_t); */
/*   if (posix_memalign((void **)&q_matrix, 64, mat_size) != 0) { */
/*     fprintf(stderr, "Memory allocation error for q_matrix\n"); */
/*     exit(-1); */
/*   } */
/*   memset(q_matrix, 0, mat_size); */
/**/
/*   int32_t *gt_arr = NULL; */
/*   int32_t ngt_arr = 0; */
/*   uint32_t check_s = 0; */
/*   clock_t START = clock(); */
/*   struct timespec STARTM, ENDM; */
/*   clock_gettime(CLOCK_MONOTONIC, &STARTM); */
/**/
/*   while (bcf_read(fp, hdr, rec) >= 0 && check_s < p.n_sites) { */
/*     bcf_unpack(rec, BCF_UN_FMT); */
/*     int nsmpl = bcf_hdr_nsamples(hdr); */
/*     int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr); */
/*     int max_ploidy = ngt / nsmpl; */
/**/
/*     for (int idx = 0; idx < nsmpl; idx++) { */
/*       int32_t *ptr = gt_arr + idx * max_ploidy; */
/*       for (int j = 0; j < max_ploidy; j++) { */
/*         if (ptr[j] == bcf_int32_vector_end) */
/*           break; */
/*         if (bcf_gt_is_missing(ptr[j])) */
/*           exit(-1); */
/*         if (bcf_gt_allele(ptr[j]) == 1) { */
/*           uint32_t q_idx = idx * 2 + j; */
/*           q_matrix[(size_t)check_s * n + q_idx] = 1; */
/*         } */
/*       } */
/*     } */
/*     check_s++; */
/*   } */
/*   clock_gettime(CLOCK_MONOTONIC, &ENDM); */
/**/
/*   double wall_time = */
/*       ENDM.tv_sec - STARTM.tv_sec + (ENDM.tv_nsec - STARTM.tv_nsec) / 1e9; */
/*   fprintf(stderr, "Queries loaded in %fs\n", wall_time); */
/**/
/*   free(gt_arr); */
/*   bcf_destroy(rec); */
/*   bcf_hdr_destroy(hdr); */
/*   bcf_close(fp); */
/**/
/*   if (check_s != p.n_sites) { */
/*     fprintf(stderr, "ERROR: PBWT sites = %u, VCF sites = %u\n", p.n_sites, */
/*             check_s); */
/*     exit(-1); */
/*   } */
/**/
/*   pbwt_col cols[3]; */
/*   memset(cols, 0, sizeof(pbwt_col) * 3); */
/*   pbwt_col *c_col = &cols[0]; */
/*   pbwt_col *next_col = &cols[1]; */
/*   pbwt_col *next_next_col = &cols[2]; */
/**/
/*   pbwt_col_deserialize(fpb, c_col); */
/*   if (p.n_sites > 1) */
/*     pbwt_col_deserialize(fpb, next_col); */
/**/
/*   q_state *states = (q_state *)malloc(n * sizeof(q_state)); */
/*   uint32_t init_c_i = p.n_haps - 1; */
/*   uint32_t init_c_r = c_col->e_pa.n - 1; */
/*   uint32_t init_c_p = c_arr_get(&c_col->e_pa, init_c_r); */
/*   uint8_t init_c_s = get_ns(c_col->zero, init_c_r); */
/**/
/*   for (size_t q = 0; q < n; q++) { */
/*     states[q].p_p = 0; */
/*     states[q].p_l = 0; */
/*     states[q].c_i = init_c_i; */
/*     states[q].c_r = init_c_r; */
/*     states[q].c_p = init_c_p; */
/*     states[q].c_s = init_c_s; */
/*   } */
/**/
/*   uint32_t i = 0; */
/**/
/* #pragma omp parallel default(none) \ */
/*     shared(n, states, q_matrix, c_col, next_col, next_next_col, p, i,
 * matches, \ */
/*                fpb) */
/*   { */
/*     while (i < p.n_sites) { */
/*       bool has_next = (i < p.n_sites - 1); */
/*       bool has_next_next = (i < p.n_sites - 2); */
/**/
/*       uint8_t *q_row = &q_matrix[(size_t)i * n]; */
/**/
/* #pragma omp single nowait */
/*       { */
/*         if (has_next_next) */
/*           pbwt_col_deserialize(fpb, next_next_col); */
/*       } */
/**/
/* #pragma omp for schedule(dynamic, 64) */
/*       for (size_t q = 0; q < n; q++) { */
/*         uint8_t q_s = q_row[q]; */
/*         q_state s = states[q]; */
/*         uint32_t cur_p, cur_l; */
/**/
/*         if (LIKELY(q_s == s.c_s)) { */
/*           cur_p = s.c_p; */
/*           cur_l = (i != 0 && s.p_p == cur_p) ? s.p_l + 1 : 1; */
/*         } else { */
/*           uint32_t t = c_arr_get(&c_col->t, s.c_r); */
/*           uint32_t p_n = c_col->p.n; */
/**/
/*           if (UNLIKELY(p_n == 1)) { */
/*             cur_p = p.n_haps; */
/*             cur_l = 0; */
/*             if (has_next) { */
/*               s.c_p = c_arr_get(&next_col->e_pa, next_col->e_pa.n - 1); */
/*               s.c_i = p.n_haps - 1; */
/*               s.c_r = next_col->e_pa.n - 1; */
/*               s.c_s = get_ns(next_col->zero, s.c_r); */
/*             } */
/*             goto check_match; */
/*           } */
/**/
/*           uint32_t d; */
/*           bool is_up = (s.c_r != 0 && ((s.c_i < t) || (s.c_r == p_n - 1)));
 */
/*           uint32_t p_i = s.c_i; */
/**/
/*           if (is_up) { */
/*             s.c_i = c_arr_get(&c_col->p, s.c_r) - 1; */
/*             s.c_p = c_arr_get(&c_col->e_pa, s.c_r - 1); */
/*             s.c_r--; */
/*             d = p_i - s.c_i; */
/*           } else { */
/*             s.c_i = c_arr_get(&c_col->p, s.c_r + 1); */
/*             s.c_p = c_arr_get(&c_col->b_pa, s.c_r + 1); */
/*             s.c_r++; */
/*             d = s.c_i - p_i; */
/*           } */
/**/
/*           if (s.p_p < p.n_haps) { */
/*             uint32_t m = */
/*                 is_up ? get_l_u(&p, s.p_p, d, i) : get_l_d(&p, s.p_p, d, i);
 */
/*             cur_l = MIN(m, s.p_l) + 1; */
/*           } else { */
/*             cur_l = 1; */
/*           } */
/*           cur_p = s.c_p; */
/*         } */
/**/
/*         if (has_next) { */
/*           s.c_i = fl(c_col, s.c_i, s.c_r); */
/*           s.c_r = get_r(next_col, s.c_i); */
/*           s.c_s = get_ns(next_col->zero, s.c_r); */
/*         } */
/**/
/*       check_match: */
/*         if (i > 0 && s.p_p != cur_p && s.p_l > 0 && s.p_l >= cur_l) { */
/*           match m = {.l = s.p_l, */
/*                      .q_n = (uint32_t)q, */
/*                      .e = i - 1, */
/*                      .h = get_haps_ms(&p, s.p_p, s.p_l, i - 1)}; */
/*           KV_PUSH_MATCH_VEC(matches[q], m); */
/*         } */
/**/
/*         if (i == p.n_sites - 1 && cur_l != 0) { */
/*           match m = {.l = cur_l, */
/*                      .q_n = (uint32_t)q, */
/*                      .e = i, */
/*                      .h = get_haps_ms(&p, cur_p, cur_l, i)}; */
/*           KV_PUSH_MATCH_VEC(matches[q], m); */
/*         } */
/**/
/*         s.p_p = cur_p; */
/*         s.p_l = cur_l; */
/*         states[q] = s; */
/*       } */
/**/
/* #pragma omp barrier */
/**/
/* #pragma omp single */
/*       { */
/*         free_pbwt_col(c_col); */
/**/
/*         pbwt_col *temp = c_col; */
/*         c_col = next_col; */
/*         next_col = next_next_col; */
/*         next_next_col = temp; */
/**/
/*         i++; */
/*       } */
/*     } */
/*   } */
/**/
/*   if (c_col) */
/*     free_pbwt_col(c_col); */
/*   if (next_col && i <= p.n_sites) */
/*     free_pbwt_col(next_col); */
/**/
/*   free(states); */
/*   free(q_matrix); */
/**/
/*   for (size_t q = 0; q < n; q++) { */
/*     match_vec_print(&matches[q]); */
/*     match_vec_free(&matches[q]); */
/*   } */
/*   free(matches); */
/**/
/*   phi_free(&p.phid); */
/*   fclose(fpb); */
/* } */
/**/
/* void pbwt_query_q(const char *pbwt_filename, const char *filename, */
/*                   int threads) { */
/*   omp_set_num_threads(threads); */
/*   clock_t START = clock(); */
/*   struct timespec STARTM, ENDM; */
/*   clock_gettime(CLOCK_MONOTONIC, &STARTM); */
/**/
/*   pbwt p; */
/*   memset(&p, 0, sizeof(pbwt)); */
/**/
/*   FILE *fpb = fopen(pbwt_filename, "rb"); */
/*   if (!fpb) { */
/*     exit(-1); */
/*   } */
/*   setvbuf(fpb, NULL, _IOFBF, 1024 * 1024 * 4); */
/**/
/*   fread(&p.n_haps, sizeof(uint32_t), 1, fpb); */
/*   fread(&p.n_sites, sizeof(uint32_t), 1, fpb); */
/*   phi_deserialize(&p.phid, fpb); */
/**/
/*   pbwt_col *pbwt_cols = (pbwt_col *)malloc(p.n_sites * sizeof(pbwt_col)); */
/*   for (uint32_t s = 0; s < p.n_sites; s++) { */
/*     pbwt_col_deserialize(fpb, &pbwt_cols[s]); */
/*   } */
/*   fclose(fpb); */
/*   clock_gettime(CLOCK_MONOTONIC, &ENDM); */
/**/
/*   double wall_time = */
/*       ENDM.tv_sec - STARTM.tv_sec + (ENDM.tv_nsec - STARTM.tv_nsec) / 1e9; */
/*   fprintf(stderr, "MuPBWT loaded in %fs\n", wall_time); */
/*   htsFile *fp = hts_open(filename, "rb"); */
/*   if (fp == NULL) { */
/*     exit(-1); */
/*   } */
/*   bcf_hdr_t *hdr = bcf_hdr_read(fp); */
/*   bcf1_t *rec = bcf_init(); */
/**/
/*   uint32_t n = bcf_hdr_nsamples(hdr) * 2; */
/**/
/*   match_vec *matches = (match_vec *)malloc(n * sizeof(match_vec)); */
/*   for (size_t q = 0; q < n; q++) { */
/*     kv_init(matches[q]); */
/*     kv_resize(match, matches[q], 256); */
/*   } */
/**/
/*   uint8_t *q_buffer[3]; */
/*   for (int k = 0; k < 3; k++) { */
/*     q_buffer[k] = (uint8_t *)malloc(n * sizeof(uint8_t)); */
/*     memset(q_buffer[k], 0, n); */
/*   } */
/**/
/*   int32_t *gt_arr = NULL; */
/*   int32_t ngt_arr = 0; */
/**/
/*   if (p.n_sites > 0) */
/*     read_bcf_row_safe(fp, hdr, rec, &gt_arr, &ngt_arr, q_buffer[0], n); */
/*   if (p.n_sites > 1) */
/*     read_bcf_row_safe(fp, hdr, rec, &gt_arr, &ngt_arr, q_buffer[1], n); */
/**/
/*   q_state *states = (q_state *)malloc(n * sizeof(q_state)); */
/*   uint32_t init_c_i = p.n_haps - 1; */
/*   uint32_t init_c_r = pbwt_cols[0].e_pa.n - 1; */
/*   uint32_t init_c_p = c_arr_get(&pbwt_cols[0].e_pa, init_c_r); */
/*   uint8_t init_c_s = get_ns(pbwt_cols[0].zero, init_c_r); */
/**/
/*   for (size_t q = 0; q < n; q++) { */
/*     states[q].p_p = 0; */
/*     states[q].p_l = 0; */
/*     states[q].c_i = init_c_i; */
/*     states[q].c_r = init_c_r; */
/*     states[q].c_p = init_c_p; */
/*     states[q].c_s = init_c_s; */
/*   } */
/**/
/*   uint32_t i = 0; */
/**/
/* #pragma omp parallel default(none) \ */
/*     shared(n, states, q_buffer, pbwt_cols, p, i, matches, fp, hdr, rec, \ */
/*                gt_arr, ngt_arr) */
/*   { */
/*     while (i < p.n_sites) { */
/*       bool has_next = (i < p.n_sites - 1); */
/*       bool has_next_next = (i < p.n_sites - 2); */
/**/
/*       uint8_t *q_row = q_buffer[i % 3]; */
/*       pbwt_col *c_col = &pbwt_cols[i]; */
/*       pbwt_col *next_col = has_next ? &pbwt_cols[i + 1] : NULL; */
/**/
/* #pragma omp single nowait */
/*       { */
/*         if (has_next_next) { */
/*           read_bcf_row_safe(fp, hdr, rec, &gt_arr, &ngt_arr, */
/*                             q_buffer[(i + 2) % 3], n); */
/*         } */
/*       } */
/**/
/* #pragma omp for schedule(dynamic, 64) */
/*       for (size_t q = 0; q < n; q++) { */
/*         if (q + 8 < n) { */
/*           __builtin_prefetch(&states[q + 8], 1, 1); */
/*           __builtin_prefetch(&q_row[q + 8], 0, 1); */
/*         } */
/**/
/*         uint8_t q_s = q_row[q]; */
/*         q_state s = states[q]; */
/*         uint32_t cur_p, cur_l; */
/**/
/*         if (LIKELY(q_s == s.c_s)) { */
/*           cur_p = s.c_p; */
/*           cur_l = (i != 0 && s.p_p == cur_p) ? s.p_l + 1 : 1; */
/*         } else { */
/*           uint32_t t = c_arr_get(&c_col->t, s.c_r); */
/*           uint32_t p_n = c_col->p.n; */
/**/
/*           if (UNLIKELY(p_n == 1)) { */
/*             cur_p = p.n_haps; */
/*             cur_l = 0; */
/*             if (has_next) { */
/*               s.c_p = c_arr_get(&next_col->e_pa, next_col->e_pa.n - 1); */
/*               s.c_i = p.n_haps - 1; */
/*               s.c_r = next_col->e_pa.n - 1; */
/*               s.c_s = get_ns(next_col->zero, s.c_r); */
/*             } */
/*             goto check_match; */
/*           } */
/**/
/*           uint32_t d; */
/*           bool is_up = (s.c_r != 0 && ((s.c_i < t) || (s.c_r == p_n - 1)));
 */
/*           uint32_t p_i = s.c_i; */
/**/
/*           if (is_up) { */
/*             s.c_i = c_arr_get(&c_col->p, s.c_r) - 1; */
/*             s.c_p = c_arr_get(&c_col->e_pa, s.c_r - 1); */
/*             s.c_r--; */
/*             d = p_i - s.c_i; */
/*           } else { */
/*             s.c_i = c_arr_get(&c_col->p, s.c_r + 1); */
/*             s.c_p = c_arr_get(&c_col->b_pa, s.c_r + 1); */
/*             s.c_r++; */
/*             d = s.c_i - p_i; */
/*           } */
/**/
/*           if (s.p_p < p.n_haps) { */
/*             uint32_t m = */
/*                 is_up ? get_l_u(&p, s.p_p, d, i) : get_l_d(&p, s.p_p, d, i);
 */
/*             cur_l = MIN(m, s.p_l) + 1; */
/*           } else { */
/*             cur_l = 1; */
/*           } */
/*           cur_p = s.c_p; */
/*         } */
/**/
/*         if (has_next) { */
/*           s.c_i = fl(c_col, s.c_i, s.c_r); */
/*           s.c_r = get_r(next_col, s.c_i); */
/*           s.c_s = get_ns(next_col->zero, s.c_r); */
/*         } */
/**/
/*       check_match: */
/*         if (i > 0 && s.p_p != cur_p && s.p_l > 0 && s.p_l >= cur_l) { */
/*           match m = {.l = s.p_l, */
/*                      .q_n = (uint32_t)q, */
/*                      .e = i - 1, */
/*                      .h = get_haps_ms(&p, s.p_p, s.p_l, i - 1)}; */
/*           KV_PUSH_MATCH_VEC(matches[q], m); */
/*         } */
/**/
/*         if (i == p.n_sites - 1 && cur_l != 0) { */
/*           match m = {.l = cur_l, */
/*                      .q_n = (uint32_t)q, */
/*                      .e = i, */
/*                      .h = get_haps_ms(&p, cur_p, cur_l, i)}; */
/*           KV_PUSH_MATCH_VEC(matches[q], m); */
/*         } */
/**/
/*         s.p_p = cur_p; */
/*         s.p_l = cur_l; */
/*         states[q] = s; */
/*       } */
/**/
/* #pragma omp barrier */
/**/
/* #pragma omp single */
/*       { */
/*         i++; */
/*       } */
/*     } */
/*   } */
/**/
/*   for (uint32_t s = 0; s < p.n_sites; s++) { */
/*     free_pbwt_col(&pbwt_cols[s]); */
/*   } */
/*   free(pbwt_cols); */
/**/
/*   free(states); */
/*   for (int k = 0; k < 3; k++) */
/*     free(q_buffer[k]); */
/*   free(gt_arr); */
/**/
/*   for (size_t q = 0; q < n; q++) { */
/*     match_vec_print(&matches[q]); */
/*     match_vec_free(&matches[q]); */
/*   } */
/*   free(matches); */
/**/
/*   bcf_destroy(rec); */
/*   bcf_hdr_destroy(hdr); */
/*   bcf_close(fp); */
/*   phi_free(&p.phid); */
/* } */
