/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

// #include "mu-pbwt/rlpbwt_int.h"
#include "mupbwt/pbwt.h"
#include <io/ref_genotype_reader.h>

#include <condition_variable>
#include <mutex>
#include <thread>

#ifdef __XSI__
#include "c_api.h"
#endif

std::map<std::string, int> mapPloidy = {{"1", 1}, {"2", 2}};

ref_genotype_reader::ref_genotype_reader(ref_haplotype_set &_H, variant_map &_V,
                                         const std::string _region,
                                         const float _sparse_maf,
                                         const bool _keep_mono)
    : H(_H), V(_V), region(_region), sparse_maf(_sparse_maf), n_ref_samples(0),
      keep_mono(_keep_mono) {
  H.sparse_maf = _sparse_maf;
}

ref_genotype_reader::~ref_genotype_reader() {}

void ref_genotype_reader::set_ploidy_ref(const int ngt_ref,
                                         const int *gt_arr_ref,
                                         const int ngt_arr_ref) {
  int max_ploidy_ref = 0;
  int n_ref_haploid = 0;

  if (gt_arr_ref == nullptr || ngt_ref == 0)
    vrb.error("Error setting ploidy while reading the first record.");

  max_ploidy_ref = ngt_ref / n_ref_samples;
  if (max_ploidy_ref < 1 || max_ploidy_ref > 2 || ngt_ref % n_ref_samples != 0)
    vrb.error(
        "Max ploidy of the reference panel cannot be set to neither 1 or 2.");

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
  H.n_ref_haps = n_ref_haploid + 2 * (n_ref_samples - n_ref_haploid);

  // Allocate some parts of haplotype set
  H.allocate();

  // Verbose
  vrb.bullet("VCF/BCF scanning [Nr=" + stb.str(n_ref_samples) +
             " (Nrh=" + stb.str(n_ref_haploid) +
             " - Nrd=" + stb.str(n_ref_samples - n_ref_haploid) + ")" + "]");
  vrb.wait("  * Reference panel parsing");
  tac.clock();
}

void ref_genotype_reader::initReader(bcf_srs_t *sr, const std::string fref,
                                     int nthreads) {
  sr->require_index = 1;
  if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1)
    vrb.error("Impossible to jump to region [" + region + "] in [" + fref +
              "]");
  if (bcf_sr_set_targets(sr, region.c_str(), 0, 0) == -1)
    vrb.error("Impossible to set target region [" + region + "] in [" + fref +
              "]");

  if (nthreads > 1)
    bcf_sr_set_threads(sr, nthreads);

  // Opening files
  if (!(bcf_sr_add_reader(sr, fref.c_str()))) {
    // we do not build an index here, as the target and reference panel could be
    // accessed in parallel
    if (sr->errnum != idx_load_failed)
      vrb.error("Failed to open file: " + fref + "");
    else
      vrb.error("Failed to load index of the file: " + fref + "");
  }
}

void ref_genotype_reader::readRefPanel(std::string fref, int nthreads,
                                       const std::string prefix_output,
                                       const std::string reg_out,
                                       bool build_mupbwt) {
  if (stb.get_extension(fref) == "xsi")
    fref += "_var.bcf";

  bcf_srs_t *sr = bcf_sr_init();
  initReader(sr, fref, nthreads);

#ifdef __XSI__
  n_ref_samples = c_xcf_nsamples(fref.c_str());
#else
  n_ref_samples = bcf_hdr_nsamples(sr->readers[0].header);
#endif

  parseRefGenotypes(sr, fref, prefix_output, reg_out, nthreads, build_mupbwt);
  bcf_sr_destroy(sr);
}

void ref_genotype_reader::parseRefGenotypes(bcf_srs_t *sr,
                                            const std::string fref,
                                            const std::string output_prefix,
                                            const std::string reg_out,
                                            int nthreads, bool build_mupbwt) {
  if (sr == nullptr)
    vrb.error("Error reading reference file");

  vrb.wait("  * Reference panel parsing");
  tac.clock();

  unsigned int i_site = 0, i_common = 0;
  int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
  bcf1_t *line_ref;
  sr->max_unpack = BCF_UN_ALL;

  unsigned int cref = 0, calt = 0;
  int line_max_ploidy = 0;
  int idx_ref_hap = 0;
  bool a = false;
  int *ptr;

  int prev_pos = -1;
  bool warning_out = false;
  const float region_span = std::max(1, V.input_stop - V.input_start);

  H.n_tot_sites = 0;
  H.n_com_sites_hq = 0;

#ifdef __XSI__
  c_xcf *c_xcf_p = c_xcf_new();
  c_xcf_add_readers(c_xcf_p, sr);
#endif

  pbwt mupbwt;
  pbwt_build_state pbwt_st;

  static const int QUEUE_DEPTH = 4;
  struct gt_slot {
    std::vector<int> buf;
    int ngt_ref = 0;
    int line_max_ploidy = 0;
  };
  std::vector<gt_slot> ring(QUEUE_DEPTH);

  std::mutex q_mutex;
  std::condition_variable q_not_empty, q_not_full;
  int q_head = 0;
  int q_tail = 0;
  int q_count = 0;
  bool q_done = false;

  std::thread pbwt_worker;

  if (build_mupbwt) {
    pbwt_build_af_init(&mupbwt, &pbwt_st, n_ref_samples);

    for (int s = 0; s < QUEUE_DEPTH; ++s)
      ring[s].buf.resize((size_t)n_ref_samples * 2);

    pbwt_worker = std::thread([&]() {
      for (;;) {
        int s;
        {
          std::unique_lock<std::mutex> lock(q_mutex);
          q_not_empty.wait(lock, [&] { return q_count > 0 || q_done; });
          if (q_count == 0)
            return;
          s = q_head;
        }
        pbwt_build_af_process_site(&mupbwt, &pbwt_st, ring[s].buf.data(),
                                   ring[s].line_max_ploidy, ploidy_ref_samples,
                                   n_ref_samples);
        {
          std::lock_guard<std::mutex> lock(q_mutex);
          q_head = (q_head + 1) % QUEUE_DEPTH;
          --q_count;
        }
        q_not_full.notify_one();
      }
    });
  }

  int rAC = 0, nAC = 0, *vAC = NULL;
  int rAN = 0, nAN = 0, *vAN = NULL;
  while (bcf_sr_next_line(sr)) {
    if (bcf_sr_get_line(sr, 0)->n_allele != 2)
      continue; // always ref

    line_ref = bcf_sr_get_line(sr, 0);

    rAC = bcf_get_info_int32(sr->readers[0].header, line_ref, "AC", &vAC, &nAC);
    rAN = bcf_get_info_int32(sr->readers[0].header, line_ref, "AN", &vAN, &nAN);
    if ((nAC != 1) || (nAN != 1))
      vrb.error(
          "VCF for reference panel needs AC/AN INFO fields to be present");
    calt = vAC[0];
    cref = (vAN[0] - vAC[0]);

    if (std::min(calt, cref) == 0 && !keep_mono) {
      if (!warning_out)
        vrb.warning(
            "Monomorphic site found [AC field] in reference panel at "
            "position: " +
            std::to_string(line_ref->pos + 1) +
            ". ALL monomorphic variants will be skipped. Please check your "
            "reference panel file. Use the --keep-monomorphic-ref-sites option "
            "to force GLIMPSE to use monomorphic sites in the reference panel. "
            "This warning is shown only once.");
      warning_out = true;
      continue;
    }

    double MAF =
        std::min(cref * 1.0f / (cref + calt), calt * 1.0f / (cref + calt));
    const bool is_common = MAF >= sparse_maf;
    H.flag_common.push_back(is_common);
    H.major_alleles.push_back(calt > cref);
    H.n_com_sites += is_common;
    H.n_rar_sites += !is_common;
    if (is_common)
      H.common2tot.push_back(H.n_tot_sites);

    int8_t line_type = bcf_get_variant_types(line_ref);
    H.n_com_sites_hq +=
        is_common && line_type == VCF_SNP && (line_ref->pos != prev_pos);
    V.push(new variant(line_ref->pos + 1, std::string(line_ref->d.id),
                       line_ref->d.allele[0], line_ref->d.allele[1], line_type,
                       V.size(), cref, calt,
                       line_type == VCF_SNP && (line_ref->pos != prev_pos)));
    prev_pos = line_ref->pos;
    H.n_tot_sites++;

#ifdef __XSI__
    ngt_ref = c_xcf_get_genotypes(c_xcf_p, 0, sr->readers[0].header, line_ref,
                                  &gt_arr_ref, &ngt_arr_ref);
#else
    ngt_ref = bcf_get_genotypes(sr->readers[0].header, line_ref, &gt_arr_ref,
                                &ngt_arr_ref);
#endif

    if (i_site == 0)
      set_ploidy_ref(ngt_ref, gt_arr_ref, ngt_arr_ref);
    line_max_ploidy = ngt_ref / n_ref_samples;
    idx_ref_hap = 0, cref = 0, calt = 0;

    if (build_mupbwt) {
      int s;
      {
        std::unique_lock<std::mutex> lock(q_mutex);
        q_not_full.wait(lock, [&] { return q_count < QUEUE_DEPTH; });
        s = q_tail;
      }
      if (ring[s].buf.size() < (size_t)ngt_ref)
        ring[s].buf.resize((size_t)ngt_ref);
      memcpy(ring[s].buf.data(), gt_arr_ref, (size_t)ngt_ref * sizeof(int));
      ring[s].ngt_ref = ngt_ref;
      ring[s].line_max_ploidy = line_max_ploidy;
      {
        std::lock_guard<std::mutex> lock(q_mutex);
        q_tail = (q_tail + 1) % QUEUE_DEPTH;
        ++q_count;
      }
      q_not_empty.notify_one();
    }

    if (is_common) {
      if (i_common >= H.HvarRef.n_rows)
        H.HvarRef.reallocate(std::max((unsigned int)(H.HvarRef.n_rows + 1),
                                      (unsigned int)(H.HvarRef.n_rows * 1.25)),
                             H.n_ref_haps);
      for (int i = 0; i < n_ref_samples; ++i) {
        ptr = gt_arr_ref + i * line_max_ploidy;
        for (int j = 0; j < ploidy_ref_samples[i]; j++) {
          a = (bcf_gt_allele(ptr[j]) == 1);
          if (a)
            H.HvarRef.set(i_common, idx_ref_hap, a);
          ++idx_ref_hap;
          a ? ++calt : ++cref;
        }
      }
    } else {
      for (int i = 0; i < n_ref_samples; ++i) {
        ptr = gt_arr_ref + i * line_max_ploidy;
        for (int j = 0; j < ploidy_ref_samples[i]; j++) {
          a = (bcf_gt_allele(ptr[j]) == 1);
          if (a != H.major_alleles[i_site])
            H.ShapRef[idx_ref_hap].push_back(i_site);
          ++idx_ref_hap;
          a ? ++calt : ++cref;
        }
      }
    }

    if ((cref != V.vec_pos[i_site]->cref) || (calt != V.vec_pos[i_site]->calt))
      vrb.error("AC/AN INFO fields in VCF are inconsistent with GT field, "
                "update the values in the VCF");
    i_common += is_common;
    i_site++;
    vrb.progress(
        "  * Reference panel parsing ",
        std::min(1.0f, (line_ref->pos + 1 - V.input_start) / region_span));
  }
  if (build_mupbwt) {
    std::lock_guard<std::mutex> lock(q_mutex);
    q_done = true;
  }
  q_not_empty.notify_one();

#ifdef __XSI__
  c_xcf_delete(c_xcf_p);
#endif
  free(vAC);
  free(vAN);
  free(gt_arr_ref);

  if (H.n_tot_sites == 0)
    vrb.error("No variants to be imputed in files");

  H.HvarRef.reallocate(H.n_com_sites, H.n_ref_haps); // trim to true final size

  vrb.bullet("Reference panel parsing done (" +
             stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");
  vrb.bullet(fref);
  vrb.bullet("VCF/BCF scanning [L=" + stb.str(H.n_tot_sites) +
             " (Lrare= " + stb.str(H.n_rar_sites) + " (" +
             stb.str(((float)H.n_rar_sites / H.n_tot_sites) * 100.0, 1) +
             "%) - Lcommon= " + stb.str(H.n_com_sites) +
             " - hq= " + stb.str(H.n_com_sites_hq) + " (" +
             stb.str(((float)H.n_com_sites_hq / H.n_tot_sites) * 100.0, 1) +
             "%))]");

  if (pbwt_worker.joinable())
    pbwt_worker.join();

  if (build_mupbwt) {
    pbwt_build_af_finalize(&mupbwt, &pbwt_st);

    auto memorize_file = output_prefix + "_" + reg_out + ".ser";
    FILE *out = fopen(memorize_file.c_str(), "wb");
    if (!out) {
      perror("fopen");
      return;
    }
    if (pbwt_serialize(out, &mupbwt) != 0) {
      fprintf(stderr, "Serialization error\n");
      fclose(out);
      return;
    }
    pbwt_print_size(&mupbwt);
    free_pbwt(&mupbwt);
  }
}
