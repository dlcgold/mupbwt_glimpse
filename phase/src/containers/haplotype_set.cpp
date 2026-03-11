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

#include "haplotype_set.h"
#include "boost/serialization/serialization.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/tmpdir.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <cmath>
#include <containers/haplotype_set.h>
#include <cstddef>
#include <iostream>
#include <math.h>
#include <objects/genotype.h>
#include <ostream>
#include <random>
#include <string>

haplotype_set::haplotype_set() {
  n_tot_sites = 0;
  n_rar_sites = 0;
  n_com_sites = 0;

  n_tot_haps = 0;
  n_tar_haps = 0;
  n_ref_haps = 0;
  n_tar_samples = 0;

  Kinit = 1000;
  K = 1000;
  Kpbwt = 2000;
  nstored = 0;
  pbwt_depth = 0;
  pbwt_modulo_cm = 0.0f;
  max_ploidy = 2;
  fploidy = 2;
  counter_sel_gf = 0;
  counter_gf = 0;
  counter_rare_restarts = 0;
}

haplotype_set::~haplotype_set() {
  n_tot_sites = 0;
  n_rar_sites = 0;
  n_com_sites = 0;

  n_tot_haps = 0;
  n_tar_haps = 0;
  n_ref_haps = 0;

  pbwt_depth = 0;
  pbwt_modulo_cm = 0;
}

void haplotype_set::allocate() {
  n_tot_haps = n_ref_haps + n_tar_haps;

  HvarRef.allocate(n_com_sites, n_ref_haps);
  ShapRef = std::vector<std::vector<int>>(n_ref_haps);
  SvarRef = std::vector<std::vector<int>>(n_tot_sites);

  HvarTar.allocate(n_com_sites, n_tar_haps);
  ShapTar = std::vector<std::vector<int>>(n_tar_haps);
  SvarTar = std::vector<std::vector<int>>(n_tot_sites);
  SindTarGL = std::vector<std::vector<int>>(n_tar_samples);

  pbwt_states = std::vector<std::vector<std::vector<int>>>(
      n_tar_samples, std::vector<std::vector<int>>());
  init_states = std::vector<std::set<int>>(n_tar_samples, std::set<int>());
  list_states = std::vector<std::vector<int>>(n_tar_haps);
}

void haplotype_set::allocate_hap_only() {
  n_tot_haps = n_ref_haps + n_tar_haps;

  SvarRef = std::vector<std::vector<int>>(n_tot_sites);
  HvarTar.allocate(n_com_sites, n_tar_haps);
  ShapTar = std::vector<std::vector<int>>(n_tar_haps);
  SvarTar = std::vector<std::vector<int>>(n_tot_sites);
  SindTarGL = std::vector<std::vector<int>>(n_tar_samples);

  pbwt_states = std::vector<std::vector<std::vector<int>>>(
      n_tar_samples, std::vector<std::vector<int>>());
  init_states = std::vector<std::set<int>>(n_tar_samples, std::set<int>());
  list_states = std::vector<std::vector<int>>(n_tar_haps);
}

void haplotype_set::initRareTar(const genotype_set &G, const variant_map &M) {
  if (Kinit == 0) {
    vrb.bullet("No init rare Tar (Kinit=0)");
    return;
  }

  tac.clock();
  stats1D statC_rare;

  std::array<float, 3> tmp;
  float sum;
  for (int i = 0; i < G.n_ind; i++) {
    const int ploidy = G.vecG[i]->ploidy;
    const int ploidyP1 = ploidy + 1;
    std::vector<unsigned char> &gls = G.vecG[i]->GL;
    std::vector<bool> &flat = G.vecG[i]->flat;

    // Clear rare alleles
    SindTarGL[i].clear();

    int counter_tot = 0;
    for (int v = 0, vc = 0; v < n_tot_sites; v++) {
      if (flat[v])
        continue;

      const int maj_gt = 2 * major_alleles[v];
      if (!flag_common[v] && !M.vec_pos[v]->LQ) {
        if (gls[ploidyP1 * v + 1] <= gls[ploidyP1 * v + maj_gt] ||
            gls[ploidyP1 * v + 2 - maj_gt] <= gls[ploidyP1 * v + maj_gt]) {
          counter_tot++;
          // let's perform calling
          tmp[0] = unphred[gls[ploidyP1 * v]];
          tmp[1] = unphred[gls[ploidyP1 * v + 1]];
          tmp[2] = ploidy > 1 ? unphred[gls[ploidyP1 * v + 2]] : 0.0f;

          sum = tmp[0] + tmp[1] + tmp[2];
          tmp[0] /= sum;
          tmp[1] /= sum;
          tmp[2] /= sum;
          // TODO //FIXME ploidy=1?
          const float af = M.vec_pos[v]->calt * 1.0 / n_ref_haps;

          tmp[0] *= (1.0f - af) * (1.0f - af);
          tmp[1] *= 2.0f * af * (1.0f - af);
          tmp[2] *= (af * af);

          sum = tmp[0] + tmp[1] + tmp[2];
          tmp[0] /= sum;
          tmp[1] /= sum;
          tmp[2] /= sum;

          if ((tmp[1] + tmp[2 - maj_gt] > tmp[maj_gt]))
            SindTarGL[i].push_back(v);
        }
      }
    }
    statC_rare.push(SindTarGL[i].size());
  }
  vrb.bullet("Init rare Tar [#non-ref gls: " + stb.str(statC_rare.mean()) +
             "] (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");
}

void haplotype_set::updateHaplotypes(const genotype_set &G) {
  tac.clock();

  for (int i = 0; i < G.n_ind; i++) {
    const int ploidy = G.vecG[i]->ploidy;
    const int hapid = tar_ind2hapid[i]; //+ n_ref_haps;

    // Clear rare alleles
    ShapTar[hapid + 0].clear();
    if (ploidy > 1)
      ShapTar[hapid + 1].clear();

    for (int v = 0, vc = 0; v < n_tot_sites; v++) {
      if (flag_common[v]) {
        HvarTar.set(vc, hapid + 0, G.vecG[i]->H0[v]);
        if (ploidy > 1)
          HvarTar.set(vc, hapid + 1, G.vecG[i]->H1[v]);
        vc++;
      } else {
        if (G.vecG[i]->H0[v] != major_alleles[v])
          ShapTar[hapid + 0].push_back(v);
        if (ploidy > 1 && G.vecG[i]->H1[v] != major_alleles[v])
          ShapTar[hapid + 1].push_back(v);
      }
    }
  }
  vrb.bullet("HAP update Tar (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) +
             "s) ");
}

void haplotype_set::updateHaplotypesTot(const genotype_set &G) {
  tac.clock();

  for (int i = 0; i < G.n_ind; i++) {
    const int ploidy = G.vecG[i]->ploidy;
    const int hapid = tar_ind2hapid[i]; //+ n_ref_haps;

    // Clear rare alleles
    ShapTar[hapid + 0].clear();
    if (ploidy > 1)
      ShapTar[hapid + 1].clear();

    for (int v = 0, vc = 0; v < n_tot_sites; v++) {

      HvarTar.set(vc, hapid + 0, G.vecG[i]->H0[v]);
      if (ploidy > 1)
        HvarTar.set(vc, hapid + 1, G.vecG[i]->H1[v]);
      vc++;
    }
  }
  vrb.bullet("HAP update Tar (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) +
             "s) ");
}

void haplotype_set::transposeRareRef() {
  tac.clock();

  // Clear previous
  for (int l = 0; l < n_tot_sites; l++)
    SvarRef[l].clear();

  // Populate
  for (int h = 0; h < n_ref_haps; h++)
    for (int r = 0; r < ShapRef[h].size(); r++)
      SvarRef[ShapRef[h][r]].push_back(h);

  vrb.bullet("RARE transpose Ref (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) +
             "s)");
}

void haplotype_set::transposeRareTar() {
  tac.clock();

  // Clear previous
  for (int l = 0; l < n_tot_sites; l++)
    SvarTar[l].clear();

  // Populate
  for (int h = 0; h < n_tar_haps; h++)
    for (int r = 0; r < ShapTar[h].size(); r++)
      SvarTar[ShapTar[h][r]].push_back(h);

  vrb.bullet("RARE transpose Tar (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) +
             "s)");
}

void haplotype_set::allocatePBWT(const int _pbwt_depth,
                                 const float _pbwt_modulo_cm,
                                 const variant_map &M, const genotype_set &G,
                                 const int _Kinit, const int _Kpbwt) {
  tac.clock();

  Kinit = _Kinit;
  Kpbwt = _Kpbwt;
  pbwt_depth = _pbwt_depth;
  pbwt_modulo_cm = _pbwt_modulo_cm;

  if (Kpbwt == 0 || Kpbwt >= n_ref_haps) {
    vrb.bullet("No PBWT allocated (Kpbwt=0 or Kpbwt >= n_ref_haps)");
    return;
  }

  pbwt_array_A = std::vector<int>(n_ref_haps);
  pbwt_array_B = std::vector<int>(n_ref_haps);
  pbwt_array_V = std::vector<int>(n_ref_haps + 1);
  pbwt_index = std::vector<int>(n_tar_haps, 0);
  f_k = std::vector<int>(n_tar_haps, 0);
  g_k = std::vector<int>(n_tar_haps, n_ref_haps);

  pbwt_small_A.reserve(n_ref_haps);
  pbwt_small_B.reserve(n_ref_haps);
  pbwt_small_V.reserve(n_ref_haps);
  pbwt_small_index = std::vector<int>(n_tar_haps);
  f_k_small = std::vector<int>(n_tar_haps);
  g_k_small = std::vector<int>(n_tar_haps);
  last_reset = std::vector<int>(n_tar_haps, 0);
  last_rare = std::vector<int>(n_tar_haps, 0);
  tar_hap = std::vector<unsigned char>(n_tar_haps);

  tar_hapid2ind = std::vector<int>(n_tar_haps);
  int idx_tar_hap = 0;
  for (int i = 0; i < n_tar_samples; ++i) {
    for (int j = 0; j < tar_ploidy[i]; ++j) {
      tar_hapid2ind[idx_tar_hap] = i;
      // cond_haps_per_hap[idx_tar_hap]=Kinit/tar_ploidy[i];
      ++idx_tar_hap;
    }
  }
  std::iota(pbwt_array_A.begin(), pbwt_array_A.end(), 0);
  pack3init();

  if (Ypacked.size() == 0)
    build_sparsePBWT(M);

  cm_pos = std::vector<float>(n_tot_sites);
  for (int i = 0; i < n_tot_sites; ++i)
    cm_pos[i] = std::max(0.0f, (float)M.vec_pos[i]->cm);

  float length = cm_pos.back() - cm_pos[0];
  if ((length / pbwt_modulo_cm) * 2 * pbwt_depth < Kpbwt) {
    while ((length / pbwt_modulo_cm) * 2 * pbwt_depth < Kpbwt &&
           pbwt_modulo_cm > 0.02)
      pbwt_modulo_cm /= 2;
    vrb.warning("Small imputation region in cM detected (for Kpbwt parameter): "
                "changing PBWT modulo parameter to: " +
                std::to_string(pbwt_modulo_cm) + " cM");
  }

  pbwt_grp.clear();
  pbwt_stored = std::vector<bool>(n_com_sites_hq, false);
  for (int l_hq = 0, l_all = 0, src = 0, tar = 0; l_all < n_com_sites;
       l_all++) {
    const int k = common2tot[l_all];
    if (M.vec_pos[k]->LQ)
      continue;
    const int tmp = (int)round(cm_pos[k] / pbwt_modulo_cm);
    if (src != tmp) {
      src = tmp;
      if (l_hq > 0)
        pbwt_grp.push_back(l_hq);
    }
    ++l_hq;
  }
  if (pbwt_grp.back() < n_com_sites_hq)
    pbwt_grp.push_back(n_com_sites_hq);
  nstored = pbwt_grp.size();

  K = pbwt_depth;
  for (int i = 0; i < n_tar_samples; ++i) {
    for (int j = 0; j < K; ++j) {
      pbwt_states[i].push_back(std::vector<int>());
      pbwt_states[i].back().reserve(tar_ploidy[i] * nstored * 2);
    }
  }
  vrb.bullet("Size PBWT (" + stb.str(Ypacked.size() / (1024 * 1024)) +
             " Mb) / pbwt-depth=" + std::to_string(K) +
             " / n_stored=" + std::to_string(nstored));
}

// void haplotype_set::matchHapsFromMuPBWT(pbwt &mupbwt, const variant_map &V,
//                                         const bool main_iteration,
//                                         std::vector<int> &uncommon_sites,
//                                         uint8_t *q_p, uint32_t n,
//                                         uint32_t n_sites, int threads) {
//   omp_set_num_threads(threads);
//   if (Kpbwt == 0 || Kpbwt >= n_ref_haps) {
//     vrb.bullet("No PBWT selection (Kpbwt=0 or Kpbwt >= n_ref_haps)");
//     return;
//   }
//
//   tac.clock();
//
//   for (int e = 0; e < n_tar_samples; e++) {
//     for (int j = 0; j < K; ++j) {
//       std::vector<int>().swap(pbwt_states[e][j]);
//     }
//   }
//
//   const size_t chunk_size = uncommon_sites.size();
//   const double thr_short = n_sites / 1000.0;
//   const double thr_medium = n_sites / 20.0;
//   uint32_t n_haps = mupbwt.n_haps;
//
// #pragma omp parallel for num_threads(threads) default(none) \
//     shared(mupbwt, q_p, n, n_sites, chunk_size, pbwt_depth, pbwt_states, \
//                tar_hapid2ind, thr_short, thr_medium, n_haps, K) \
//     schedule(dynamic)
//   for (size_t q = 0; q < n; ++q) {
//     uint32_t target_ind = tar_hapid2ind[q];
//     std::vector<double> haplo_score(n_haps, 0.0);
//     const uint8_t *query = &q_p[q * n_sites];
//
//     uint32_t p = 0, p_p = 0;
//     uint32_t l = 0, p_l = 0;
//
//     pbwt_col *all_cols = mupbwt.cols.a;
//     pbwt_col *c_col = &all_cols[0];
//     uint32_t c_p = c_arr_get(&c_col->e_pa, c_col->e_pa.n - 1);
//     uint32_t c_i = n_haps - 1;
//     uint32_t c_r = c_col->e_pa.n - 1;
//     uint8_t c_s = get_ns(c_col->zero, c_r);
//     uint32_t p_i = c_i;
//
//     for (size_t i = 0; i < n_sites; i++) {
//       c_col = &all_cols[i];
//       pbwt_col *next_col = (i < n_sites - 1) ? &all_cols[i + 1] : NULL;
//       uint8_t q_s = query[i];
//
//       if (LIKELY(q_s == c_s)) {
//         p = c_p;
//         l = (i != 0 && p_p == p) ? p_l + 1 : 1;
//       } else {
//         uint32_t t = c_arr_get(&c_col->t, c_r);
//         uint32_t p_n = c_col->p.n;
//
//         if (UNLIKELY(p_n == 1)) {
//           p = n_haps;
//           l = 0;
//           if (next_col) {
//             c_p = c_arr_get(&next_col->e_pa, next_col->e_pa.n - 1);
//             c_i = n_haps - 1;
//             c_r = next_col->e_pa.n - 1;
//             c_s = get_ns(next_col->zero, c_r);
//           }
//           goto check_match;
//         }
//
//         uint32_t d;
//         bool is_up = (c_r != 0 && ((c_i < t) || (c_r == p_n - 1)));
//
//         if (is_up) {
//           p_i = c_i;
//           c_i = c_arr_get(&c_col->p, c_r) - 1;
//           c_p = c_arr_get(&c_col->e_pa, c_r - 1);
//           c_r--;
//           d = p_i - c_i;
//         } else {
//           p_i = c_i;
//           c_i = c_arr_get(&c_col->p, c_r + 1);
//           c_p = c_arr_get(&c_col->b_pa, c_r + 1);
//           c_r++;
//           d = c_i - p_i;
//         }
//
//         if (d < 90) {
//           uint32_t m =
//               is_up ? get_l_u(&mupbwt, p_p, d, i) : get_l_d(&mupbwt, p_p, d,
//               i);
//           l = std::min((uint32_t)m, p_l) + 1;
//         } else {
//           l = get_l(&mupbwt, query, i, c_i);
//         }
//         p = c_p;
//       }
//
//       if (next_col) {
//         c_i = fl(c_col, c_i, c_r);
//         c_r = get_r(next_col, c_i);
//         c_s = get_ns(next_col->zero, c_r);
//       }
//
//     check_match:
//       if (i > 0 && p_p != p && p_l > 0 && p_l >= l) {
//         int_vec h_vec;
//         double contrib = static_cast<double>(p_l) / n_sites;
//         double weight;
//
//         if (p_l < thr_short) {
//           h_vec = get_haps_fixed_depth(&mupbwt, p_p, i - 1, K / 2);
//           weight = contrib;
//         } else if (p_l < thr_medium) {
//           h_vec = get_haps_fixed_depth(&mupbwt, p_p, i - 1, K);
//           weight = contrib;
//         } else {
//           h_vec = get_haps_ms(&mupbwt, p_p, p_l, i - 1);
//           weight = contrib * contrib;
//         }
//
//         for (size_t k = 0; k < h_vec.n; k++) {
//           haplo_score[h_vec.a[k]] += weight;
//         }
//         free(h_vec.a);
//       }
//
//       if (i == n_sites - 1 && l != 0) {
//         int_vec h_vec;
//         double contrib = static_cast<double>(l) / n_sites;
//         double weight;
//
//         if (l < thr_short) {
//           h_vec = get_haps_fixed_depth(&mupbwt, p, i, K / 2);
//           weight = contrib;
//         } else if (l < thr_medium) {
//           h_vec = get_haps_fixed_depth(&mupbwt, p, i, K);
//           weight = contrib;
//         } else {
//           h_vec = get_haps_ms(&mupbwt, p, l, i);
//           weight = contrib * contrib;
//         }
//
//         for (size_t k = 0; k < h_vec.n; k++) {
//           haplo_score[h_vec.a[k]] += weight;
//         }
//         free(h_vec.a);
//       }
//
//       p_p = p;
//       p_l = l;
//     }
//
//     std::vector<std::pair<int, double>> haplo_vec;
//     haplo_vec.reserve(n_haps);
//     double tot_s = 0.0;
//
//     for (uint32_t h = 0; h < n_haps; h++) {
//       if (haplo_score[h] > 0.0) {
//         double final_score = (haplo_score[h] > 1.0) ? 1.0 : haplo_score[h];
//         haplo_vec.emplace_back(h, final_score);
//         tot_s += final_score;
//       }
//     }
//
//     double avg_s = haplo_vec.empty() ? 0.0 : tot_s / haplo_vec.size();
//
//     std::sort(haplo_vec.begin(), haplo_vec.end(),
//               [](const auto &a, const auto &b) { return a.second > b.second;
//               });
//
//     size_t idx = 0;
//     int c_p_count = 0;
// #pragma omp critical
//     {
//       for (size_t d = 0; d < pbwt_depth && idx < haplo_vec.size(); d++) {
//         if (haplo_vec[idx].second < avg_s)
//           break;
//
//         size_t end = std::min(idx + chunk_size, haplo_vec.size());
//         for (size_t i = idx; i < end; i++) {
//           pbwt_states[target_ind][d].push_back(haplo_vec[i].first);
//           c_p_count++;
//         }
//         idx = end;
//       }
//
//       if (c_p_count == 0 && !haplo_vec.empty()) {
//
//         double check = haplo_vec[0].second;
//         size_t cc = 0;
//         for (const auto &[h, s] : haplo_vec) {
//           pbwt_states[target_ind][pbwt_depth - 1].push_back(h);
//           if (cc >= chunk_size || (check - s) > 0.1)
//             break;
//           check = s;
//           cc++;
//         }
//       }
//     }
//   }
//
//   vrb.bullet("Mu-PBWT selection (" + stb.str(tac.rel_time() * 1.0 / 1000, 2)
//   +
//              "s)");
// }
// void haplotype_set::matchHapsFromMuPBWT(pbwt &mupbwt, const variant_map &V,
//                                         const bool main_iteration,
//                                         std::vector<int> &uncommon_sites,
//                                         uint8_t *q_p, uint32_t n,
//                                         uint32_t n_sites, int threads) {
//   omp_set_num_threads(threads);
//
//   if (Kpbwt == 0 || Kpbwt >= n_ref_haps) {
//     vrb.bullet("No PBWT selection (Kpbwt=0 or Kpbwt >= n_ref_haps)");
//     return;
//   }
//
//   tac.clock();
//
//   for (int e = 0; e < n_tar_samples; e++) {
//     for (int j = 0; j < K; ++j) {
//       std::vector<int>().swap(pbwt_states[e][j]);
//     }
//   }
//
//   const size_t chunk_size = uncommon_sites.size();
//   const double thr_short = n_sites / 1000.0;
//   const double thr_medium = n_sites / 20.0;
//   uint32_t n_haps = mupbwt.n_haps;
//
// #pragma omp parallel for num_threads(threads) default(none) \
//     shared(mupbwt, q_p, n, n_sites, chunk_size, pbwt_depth, pbwt_states, \
//                tar_hapid2ind, thr_short, thr_medium, n_haps, K) \
//     schedule(dynamic)
//   for (size_t q = 0; q < n; ++q) {
//     uint32_t target_ind = tar_hapid2ind[q];
//     std::vector<double> haplo_score(n_haps, 0.0);
//     std::vector<uint32_t> active_haplos;
//     active_haplos.reserve(1000);
//
//     const uint8_t *query = &q_p[q * n_sites];
//
//     uint32_t p = 0, p_p = 0;
//     uint32_t l = 0, p_l = 0;
//
//     pbwt_col *all_cols = mupbwt.cols.a;
//     pbwt_col *c_col = &all_cols[0];
//     uint32_t c_p = c_arr_get(&c_col->e_pa, c_col->e_pa.n - 1);
//     uint32_t c_i = n_haps - 1;
//     uint32_t c_r = c_col->e_pa.n - 1;
//     uint8_t c_s = get_ns(c_col->zero, c_r);
//     uint32_t p_i = c_i;
//
//     for (size_t i = 0; i < n_sites; i++) {
//       c_col = &all_cols[i];
//       pbwt_col *next_col = (i < n_sites - 1) ? &all_cols[i + 1] : NULL;
//       uint8_t q_s = query[i];
//
//       if (LIKELY(q_s == c_s)) {
//         p = c_p;
//         l = (i != 0 && p_p == p) ? p_l + 1 : 1;
//       } else {
//         uint32_t t = c_arr_get(&c_col->t, c_r);
//         uint32_t p_n = c_col->p.n;
//
//         if (UNLIKELY(p_n == 1)) {
//           p = n_haps;
//           l = 0;
//           if (next_col) {
//             c_p = c_arr_get(&next_col->e_pa, next_col->e_pa.n - 1);
//             c_i = n_haps - 1;
//             c_r = next_col->e_pa.n - 1;
//             c_s = get_ns(next_col->zero, c_r);
//           }
//           goto check_match;
//         }
//
//         uint32_t d;
//         bool is_up = (c_r != 0 && ((c_i < t) || (c_r == p_n - 1)));
//
//         if (is_up) {
//           p_i = c_i;
//           c_i = c_arr_get(&c_col->p, c_r) - 1;
//           c_p = c_arr_get(&c_col->e_pa, c_r - 1);
//           c_r--;
//           d = p_i - c_i;
//         } else {
//           p_i = c_i;
//           c_i = c_arr_get(&c_col->p, c_r + 1);
//           c_p = c_arr_get(&c_col->b_pa, c_r + 1);
//           c_r++;
//           d = c_i - p_i;
//         }
//
//         if (d < 90) {
//           uint32_t m =
//               is_up ? get_l_u(&mupbwt, p_p, d, i) : get_l_d(&mupbwt, p_p, d,
//               i);
//           l = std::min((uint32_t)m, p_l) + 1;
//         } else {
//           l = get_l(&mupbwt, query, i, c_i);
//         }
//         p = c_p;
//       }
//
//       if (next_col) {
//         c_i = fl(c_col, c_i, c_r);
//         c_r = get_r(next_col, c_i);
//         c_s = get_ns(next_col->zero, c_r);
//       }
//
//     check_match:
//       if (i > 0 && p_p != p && p_l > 0 && p_l >= l) {
//         int_vec h_vec;
//         double contrib = static_cast<double>(p_l) / n_sites;
//         double weight;
//
//         if (p_l < thr_short) {
//           h_vec = get_haps_fixed_depth(&mupbwt, p_p, i - 1, K / 2);
//           weight = contrib;
//         } else if (p_l < thr_medium) {
//           h_vec = get_haps_fixed_depth(&mupbwt, p_p, i - 1, K);
//           weight = contrib;
//         } else {
//           h_vec = get_haps_ms(&mupbwt, p_p, p_l, i - 1);
//           weight = contrib * contrib;
//         }
//
//         for (size_t k = 0; k < h_vec.n; k++) {
//           uint32_t h_id = h_vec.a[k];
//           if (haplo_score[h_id] == 0.0) {
//             active_haplos.push_back(h_id);
//           }
//           haplo_score[h_id] += weight;
//         }
//         free(h_vec.a);
//       }
//
//       if (i == n_sites - 1 && l != 0) {
//         int_vec h_vec;
//         double contrib = static_cast<double>(l) / n_sites;
//         double weight;
//
//         if (l < thr_short) {
//           h_vec = get_haps_fixed_depth(&mupbwt, p, i, K / 2);
//           weight = contrib;
//         } else if (l < thr_medium) {
//           h_vec = get_haps_fixed_depth(&mupbwt, p, i, K);
//           weight = contrib;
//         } else {
//           h_vec = get_haps_ms(&mupbwt, p, l, i);
//           weight = contrib * contrib;
//         }
//
//         for (size_t k = 0; k < h_vec.n; k++) {
//           uint32_t h_id = h_vec.a[k];
//           if (haplo_score[h_id] == 0.0) {
//             active_haplos.push_back(h_id);
//           }
//           haplo_score[h_id] += weight;
//         }
//         free(h_vec.a);
//       }
//
//       p_p = p;
//       p_l = l;
//     }
//
//     std::vector<std::pair<int, double>> haplo_vec;
//     haplo_vec.reserve(active_haplos.size());
//     double tot_s = 0.0;
//
//     for (uint32_t h : active_haplos) {
//       double final_score = std::min(haplo_score[h], 1.0);
//       haplo_vec.emplace_back(h, final_score);
//       tot_s += final_score;
//     }
//
//     double avg_s = haplo_vec.empty() ? 0.0 : tot_s / haplo_vec.size();
//
//     size_t elements_needed = pbwt_depth * chunk_size;
//     size_t top_k = std::min(haplo_vec.size(), elements_needed);
//
//     if (top_k > 0) {
//       std::partial_sort(
//           haplo_vec.begin(), haplo_vec.begin() + top_k, haplo_vec.end(),
//           [](const auto &a, const auto &b) { return a.second > b.second; });
//     }
//
//     std::vector<std::vector<int>> local_states(pbwt_depth);
//     size_t idx = 0;
//     int c_p_count = 0;
//
//     for (size_t d = 0; d < pbwt_depth && idx < top_k; d++) {
//       if (haplo_vec[idx].second < avg_s)
//         break;
//
//       size_t end = std::min(idx + chunk_size, top_k);
//       for (size_t i = idx; i < end; i++) {
//         local_states[d].push_back(haplo_vec[i].first);
//         c_p_count++;
//       }
//       idx = end;
//     }
//
//     if (c_p_count == 0 && !haplo_vec.empty()) {
//       double check = haplo_vec[0].second;
//       size_t cc = 0;
//       for (const auto &[h, s] : haplo_vec) {
//         local_states[pbwt_depth - 1].push_back(h);
//         if (cc >= chunk_size || (check - s) > 0.1)
//           break;
//         check = s;
//         cc++;
//       }
//     }
//
// #pragma omp critical
//     {
//       for (size_t d = 0; d < pbwt_depth; d++) {
//         if (!local_states[d].empty()) {
//           pbwt_states[target_ind][d].insert(pbwt_states[target_ind][d].end(),
//                                             local_states[d].begin(),
//                                             local_states[d].end());
//         }
//       }
//     }
//   }
//
//   vrb.bullet("Mu-PBWT selection (" + stb.str(tac.rel_time() * 1.0 / 1000, 2)
//   +
//              "s)");
// }

void haplotype_set::matchHapsFromMuPBWT(pbwt &mupbwt, const variant_map &V,
                                        const bool main_iteration,
                                        std::vector<int> &uncommon_sites,
                                        uint8_t *q_p, uint32_t n,
                                        uint32_t n_sites, int threads) {
  omp_set_num_threads(threads);

  if (Kpbwt == 0 || Kpbwt >= n_ref_haps) {
    vrb.bullet("No PBWT selection (Kpbwt=0 or Kpbwt >= n_ref_haps)");
    return;
  }

  tac.clock();

  for (int e = 0; e < n_tar_samples; e++) {
    for (int j = 0; j < K; ++j) {
      std::vector<int>().swap(pbwt_states[e][j]);
    }
  }

  const size_t chunk_size = uncommon_sites.size();
  const double thr_short = n_sites / 1000.0;
  const double thr_medium = n_sites / 20.0;
  uint32_t n_haps = mupbwt.n_haps;

#pragma omp parallel num_threads(threads) default(none)                        \
    shared(mupbwt, q_p, n, n_sites, chunk_size, pbwt_depth, pbwt_states,       \
               tar_hapid2ind, thr_short, thr_medium, n_haps, K)
  {
    std::vector<double> haplo_score(n_haps, 0.0);
    std::vector<uint32_t> active_haplos;
    active_haplos.reserve(2000);
    std::vector<std::vector<int>> local_states(pbwt_depth);

#pragma omp for schedule(dynamic)
    for (size_t q = 0; q < n; ++q) {
      uint32_t target_ind = tar_hapid2ind[q];

      for (uint32_t h : active_haplos) {
        haplo_score[h] = 0.0;
      }
      active_haplos.clear();
      for (auto &vec : local_states) {
        vec.clear();
      }

      const uint8_t *query = &q_p[q * n_sites];

      uint32_t p = 0, p_p = 0;
      uint32_t l = 0, p_l = 0;

      pbwt_col *all_cols = mupbwt.cols.a;
      pbwt_col *c_col = &all_cols[0];
      uint32_t c_p = c_arr_get(&c_col->e_pa, c_col->e_pa.n - 1);
      uint32_t c_i = n_haps - 1;
      uint32_t c_r = c_col->e_pa.n - 1;
      uint8_t c_s = get_ns(c_col->zero, c_r);
      uint32_t p_i = c_i;

      auto apply_score_inline = [&](uint32_t start_p, uint32_t match_len,
                                    uint32_t col) {
        double contrib = static_cast<double>(match_len) / n_sites;
        double weight = contrib;
        uint32_t depth_limit = 0;
        bool exact_only = false;

        if (match_len < thr_short) {
          depth_limit = K / 2;
        } else if (match_len < thr_medium) {
          depth_limit = K;
        } else {
          weight = contrib * contrib;
          exact_only = true;
        }

        auto add_score = [&](uint32_t h_id) {
          if (haplo_score[h_id] == 0.0) {
            active_haplos.push_back(h_id);
          }
          haplo_score[h_id] += weight;
        };

        add_score(start_p);

        uint32_t t_p = start_p;
        uint32_t steps = 0;
        while (true) {
          uint32_t d_p = phi_inv_f(&mupbwt.phid, t_p, col + 1);
          if (d_p == n_haps)
            break;

          if (exact_only) {
            if (phi_l(&mupbwt.phid, d_p, col + 1) >= match_len) {
              add_score(d_p);
              t_p = d_p;
            } else {
              break;
            }
          } else {
            if (steps < depth_limit) {
              add_score(d_p);
              t_p = d_p;
              steps++;
            } else {
              break;
            }
          }
        }

        t_p = start_p;
        steps = 0;
        while (true) {
          uint32_t u_p = phi_f(&mupbwt.phid, t_p, col + 1);
          if (u_p == n_haps)
            break;

          if (exact_only) {
            if (phi_l(&mupbwt.phid, t_p, col + 1) >= match_len) {
              add_score(u_p);
              t_p = u_p;
            } else {
              break;
            }
          } else {
            if (steps < depth_limit) {
              add_score(u_p);
              t_p = u_p;
              steps++;
            } else {
              break;
            }
          }
        }
      };

      for (size_t i = 0; i < n_sites; i++) {
        c_col = &all_cols[i];
        pbwt_col *next_col = (i < n_sites - 1) ? &all_cols[i + 1] : NULL;
        uint8_t q_s = query[i];

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
            uint32_t m = is_up ? get_l_u(&mupbwt, p_p, d, i)
                               : get_l_d(&mupbwt, p_p, d, i);
            l = std::min((uint32_t)m, p_l) + 1;
          } else {
            l = get_l(&mupbwt, query, i, c_i);
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
          apply_score_inline(p_p, p_l, i - 1);
        }

        if (i == n_sites - 1 && l != 0) {
          apply_score_inline(p, l, i);
        }

        p_p = p;
        p_l = l;
      }

      std::vector<std::pair<int, double>> haplo_vec;
      haplo_vec.reserve(active_haplos.size());
      double tot_s = 0.0;

      for (uint32_t h : active_haplos) {
        double final_score = std::min(haplo_score[h], 1.0);
        haplo_vec.emplace_back(h, final_score);
        tot_s += final_score;
      }

      double avg_s = haplo_vec.empty() ? 0.0 : tot_s / haplo_vec.size();
      size_t elements_needed = pbwt_depth * chunk_size;
      size_t top_k = std::min(haplo_vec.size(), elements_needed);

      if (top_k > 0) {
        std::partial_sort(
            haplo_vec.begin(), haplo_vec.begin() + top_k, haplo_vec.end(),
            [](const auto &a, const auto &b) { return a.second > b.second; });
      }

      if (top_k > 0) {
        size_t idx = 0;
        for (size_t d = 0; d < pbwt_depth && idx < top_k; d++) {
          if (haplo_vec[idx].second < avg_s)
            break;
          size_t end = std::min(idx + chunk_size, top_k);
          for (size_t i = idx; i < end; i++) {
            local_states[d].push_back(haplo_vec[i].first);
          }
          idx = end;
        }
      } else {
        uint32_t t_p = p;
        uint32_t steps = 0;
        local_states[pbwt_depth - 1].push_back(p);

        while (steps < chunk_size) {
          uint32_t d_p = phi_inv_f(&mupbwt.phid, t_p, n_sites);
          if (d_p == n_haps)
            break;
          local_states[pbwt_depth - 1].push_back(d_p);
          t_p = d_p;
          steps++;
        }
      }

#pragma omp critical
      {
        for (size_t d = 0; d < pbwt_depth; d++) {
          if (!local_states[d].empty()) {
            pbwt_states[target_ind][d].insert(pbwt_states[target_ind][d].end(),
                                              local_states[d].begin(),
                                              local_states[d].end());
          }
        }
      }
    }
  }

  vrb.bullet("Mu-PBWT selection (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) +
             "s)");
}

std::vector<unsigned int>
haplotype_set::extract_random(const std::vector<unsigned int> &input,
                              std::size_t n) {
  std::vector<unsigned int> r;
  auto elems = input;
  if (elems.empty() || n == 0)
    return r;

  std::random_device rd;
  std::mt19937 gen(rd());

  if (n >= elems.size()) {
    return elems;
  }

  std::shuffle(elems.begin(), elems.end(), gen);
  for (std::size_t i = 0; i < n; ++i) {
    r.push_back(elems[i]);
  }

  return r;
}

void haplotype_set::matchHapsFromCompressedPBWTSmall(
    const variant_map &M, const bool main_iteration) {
  if (Kpbwt == 0 || Kpbwt >= n_ref_haps) {
    vrb.bullet("No PBWT selection (Kpbwt=0 or Kpbwt >= n_ref_haps)");
    return;
  }

  std::cout << std::endl;
  tac.clock();
  const unsigned char *pY = &Ypacked[0];

  std::fill(last_reset.begin(), last_reset.end(), 0);
  std::fill(last_rare.begin(), last_rare.end(), 0);
  std::iota(pbwt_array_A.begin(), pbwt_array_A.end(), 0);

  stats1D size_comp;
  stats1D nrare;

  std::fill(pbwt_stored.begin(), pbwt_stored.end(), false);
  for (int idx = 0, loffset = 0; idx < pbwt_grp.size(); ++idx) {
    const int state = rng.getInt(loffset, pbwt_grp[idx] - 1);
    pbwt_stored[state] = true;
    loffset = pbwt_grp[idx];
  }

  for (int e = 0; e < n_tar_haps; e++) {
    f_k[e] = 0;
    g_k[e] = n_ref_haps;
    pbwt_index[e] = rng.getInt(0, n_ref_haps - 1);
  }
  for (int e = 0; e < n_tar_samples; e++)
    for (int j = 0; j < K; ++j)
      pbwt_states[e][j].clear();

  int ref_rac_l_com = 0, prev_ref_rac_l_com = 0;
  int ref_rac_l_rare = 0, prev_ref_rac_l_rare = 0;

  length_sel_mod.clear();
  length_sel_gf.clear();
  length_sel_gf_rare.clear();
  counter_gf = 0;
  counter_sel_gf = 0;
  counter_rare_restarts = 0;
  int l_hq = 0;
  int last_k = -1;
  int l_all = 0;
  for (int k = 0; k < n_tot_sites; k++) {
    if (M.vec_pos[k]->LQ) {
      l_all += flag_common[k];
      continue;
    }

    if (flag_common[k]) {
      if (A_small_idx[l_hq].size() > 0 && !(last_k >= 0 && flag_common[last_k]))
        init_common(k, l_hq, prev_ref_rac_l_com);
      ref_rac_l_com = M.vec_pos[k]->cref;
      read_full_pbwt_av(pY, ref_rac_l_com);
      select_common_pd_fg(k, l_hq, l_all, ref_rac_l_com, prev_ref_rac_l_com);
      prev_ref_rac_l_com = ref_rac_l_com;
      ++l_all;
      ++l_hq;
    } else {
      if (last_k < 0 || flag_common[last_k])
        init_rare(M, k, l_hq);
      ref_rac_l_rare = major_alleles[k]
                           ? M.vec_pos[k]->cref
                           : A_small_idx[l_hq].size() - M.vec_pos[k]->calt;
      read_small_pbwt_av(pY, ref_rac_l_rare, rareTarHaps.size());
      select_rare_pd_fg(k, ref_rac_l_rare);
      prev_ref_rac_l_rare = ref_rac_l_rare;
    }
    last_k = k;
  }

  vrb.bullet("sparsePBWT selection (" +
             stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");
  // vrb.bullet("Nrare: " + to_string(nrare.mean()) + " / Rare restarts: " +
  // to_string(counter_rare_restarts/(1.0*n_tar_samples)));
  // vrb.bullet(to_string(length_sel_mod.mean()) + " cM / count_sel_gf: " +
  // to_string(counter_sel_gf/(1.0*n_tar_samples)) + " / " +
  // to_string(length_sel_gf.mean()) + " cM / count_sel_gf_rare: " +
  // to_string(length_sel_gf_rare.mean()));
}

void haplotype_set::read_full_pbwt_av(const unsigned char *&pY,
                                      const int ref_rac_l) {
  unsigned char z;
  int u = 0;
  int v = ref_rac_l;
  int m = 0;
  int n;
  int mm;
  const int size_a = n_ref_haps;
  const std::array<int *, 2> occ = {&u, &v};
  pbwt_array_V[0] = 0;
  while (m < size_a) {
    z = *pY++;              // read the encoded symbol
    n = p3decode[z & 0x7f]; // get how many chars are encoded
    mm = m + n;             // get limit of this block
    z >>= 7;                // get the symbol itself
    std::copy(pbwt_array_A.begin() + m, pbwt_array_A.begin() + mm,
              pbwt_array_B.begin() + *occ[z]); // update B, (next A)
    z ? std::iota(pbwt_array_V.begin() + m + 1, pbwt_array_V.begin() + mm + 1,
                  v - ref_rac_l + 1)
      : std::fill(
            pbwt_array_V.begin() + m + 1, pbwt_array_V.begin() + mm + 1,
            v - ref_rac_l); // we pick V as we expect more fill() than iota
    m = mm;                 // update position in PPA for next loop
    *occ[z] += n;           // update FM
  }
  pbwt_array_B.swap(pbwt_array_A); // switch A and B
}

void haplotype_set::read_small_pbwt_av(const unsigned char *&pY,
                                       const int ref_rac_l,
                                       const bool update_v) {
  unsigned char z;
  int u = 0;
  int v = ref_rac_l;
  int m = 0;
  int n;
  int mm;
  const int size_a = pbwt_small_A.size();
  const std::array<int *, 2> occ = {&u, &v};
  pbwt_small_V[0] = 0;
  while (m < size_a) {
    z = *pY++;              // read the encoded symbol
    n = p3decode[z & 0x7f]; // get how many chars are encoded
    z >>= 7;                // get the symbol itself
    mm = m + n;             // get limit of this block
    std::copy(pbwt_small_A.begin() + m, pbwt_small_A.begin() + mm,
              pbwt_small_B.begin() + *occ[z]); // update B, (next A)
    if (update_v)
      z ? std::iota(pbwt_small_V.begin() + m + 1, pbwt_small_V.begin() + mm + 1,
                    v - ref_rac_l + 1)
        : std::fill(pbwt_small_V.begin() + m + 1, pbwt_small_V.begin() + mm + 1,
                    v - ref_rac_l);
    m = mm;       // update position in PPA for next loop
    *occ[z] += n; // update FM
  }
  pbwt_small_B.swap(pbwt_small_A); // switch A and B
}

void haplotype_set::select_common_pd_fg(const int k, const int l_hq,
                                        const int l_all, const int ref_rac_l,
                                        const int prev_ref_rac_l) {
  for (int htr = 0; htr < n_tar_haps; htr++) {
    bool reset = false;
    int f_dash, g_dash, idx;
    const unsigned char prev_hap = tar_hap[htr];
    tar_hap[htr] = HvarTar.get(l_all, htr);

    idx =
        tar_hap[htr] * (ref_rac_l + pbwt_array_V[pbwt_index[htr]]) +
        (1 - tar_hap[htr]) * (pbwt_index[htr] - pbwt_array_V[pbwt_index[htr]]);
    f_dash = tar_hap[htr] * (ref_rac_l + pbwt_array_V[f_k[htr]]) +
             (1 - tar_hap[htr]) * (f_k[htr] - pbwt_array_V[f_k[htr]]);
    g_dash = tar_hap[htr] * (ref_rac_l + pbwt_array_V[g_k[htr]]) +
             (1 - tar_hap[htr]) * (g_k[htr] - pbwt_array_V[g_k[htr]]);

    if (g_dash <= f_dash) {
      ++counter_gf;
      if (pbwt_stored[l_hq] || cm_pos[k - 1] - cm_pos[last_reset[htr]] >
                                   pbwt_modulo_cm / 2.0f) // only long matches
      {
        counter_sel_gf++;
        selectK(htr, k - 1, prev_ref_rac_l, pbwt_array_B, K, prev_hap);
        length_sel_gf.push(cm_pos[k - 1] - cm_pos[last_reset[htr]]);
      }
      f_dash = tar_hap[htr] * ref_rac_l;
      g_dash = tar_hap[htr] * n_ref_haps + (1 - tar_hap[htr]) * ref_rac_l;
      // f_dash=0;
      // g_dash=n_ref_haps;
      last_reset[htr] = k;
      reset = true;
    }

    pbwt_index[htr] = idx;
    f_k[htr] = f_dash;
    g_k[htr] = g_dash;

    if (!reset && pbwt_stored[l_hq]) {
      selectK(htr, k, ref_rac_l, pbwt_array_A, K, tar_hap[htr]);
      length_sel_mod.push(cm_pos[k] - cm_pos[last_reset[htr]]);
      //			f_k[htr]=0;
      //			g_k[htr]=n_ref_haps;
      f_k[htr] = tar_hap[htr] * ref_rac_l;
      g_k[htr] = tar_hap[htr] * n_ref_haps + (1 - tar_hap[htr]) * ref_rac_l;
      last_reset[htr] = k;
    }
  }
}

void haplotype_set::select_rare_pd_fg(const int k, const int ref_rac_l) {
  if (rareTarHaps.size() == 0)
    return;

  const int n_haps_small = pbwt_small_A.size();
  const int n_tar_rare = rareTarHaps.size();

  int id_rare = 0;
  int id_minor = 0;
  for (int j = 0; j < SvarTar[k].size(); ++j) {
    const int htr = SvarTar[k][j];
    last_rare[htr] = k;
    const int small_idx = ref_rac_l + pbwt_small_V[pbwt_small_index[htr]];
    int f_dash = (ref_rac_l + pbwt_small_V[f_k_small[htr]]);
    int g_dash = (ref_rac_l + pbwt_small_V[g_k_small[htr]]);
    if (g_dash <= f_dash) {
      f_dash = ref_rac_l;
      g_dash = n_haps_small;
      last_reset[htr] = k;
      ++counter_rare_restarts;
    }
    pbwt_small_index[htr] = small_idx;
    f_k_small[htr] = f_dash;
    g_k_small[htr] = g_dash;

    while (rareTarHaps[id_rare] < htr) {
      const int htr2 = rareTarHaps[id_rare];
      const int small_idx2 =
          pbwt_small_index[htr2] - pbwt_small_V[pbwt_small_index[htr2]];
      int f_dash2 = (f_k_small[htr2] - pbwt_small_V[f_k_small[htr2]]);
      int g_dash2 = (g_k_small[htr2] - pbwt_small_V[g_k_small[htr2]]);

      if (g_dash2 <= f_dash2) {
        f_dash2 = 0;
        g_dash2 = ref_rac_l;
        last_reset[htr2] = k;
        ++counter_rare_restarts;
      }
      pbwt_small_index[htr2] = small_idx2;
      f_k_small[htr2] = f_dash2;
      g_k_small[htr2] = g_dash2;
      ++id_rare;
    }
    ++id_rare;
  }

  while (id_rare < n_tar_rare) {
    const int htr2 = rareTarHaps[id_rare];
    const int small_idx2 =
        pbwt_small_index[htr2] - pbwt_small_V[pbwt_small_index[htr2]];
    int f_dash2 = (f_k_small[htr2] - pbwt_small_V[f_k_small[htr2]]);
    int g_dash2 = (g_k_small[htr2] - pbwt_small_V[g_k_small[htr2]]);

    if (g_dash2 <= f_dash2) {
      f_dash2 = 0;
      g_dash2 = ref_rac_l;
      last_reset[htr2] = k;
      ++counter_rare_restarts;
    }
    pbwt_small_index[htr2] = small_idx2;
    f_k_small[htr2] = f_dash2;
    g_k_small[htr2] = g_dash2;
    ++id_rare;
  }
}

void haplotype_set::init_common(const int k, const int l,
                                const int prev_ref_rac_l_com) {
  if (A_small_idx[l].size() <= 0)
    return;

  int j = 0;
  for (int htr = 0; htr < n_tar_haps; ++htr) {
    int f_dash, g_dash, idx;
    if (j < rareTarHaps.size() &&
        htr == rareTarHaps[j]) // rare2common //TODO check this
    {
      idx =
          n_ref_haps - A_small_idx[l].size() + pbwt_small_index[rareTarHaps[j]];
      f_dash =
          last_reset[htr] >= last_rare[htr]
              ? 0
              : n_ref_haps - A_small_idx[l].size() + f_k_small[rareTarHaps[j]];
      g_dash = n_ref_haps - A_small_idx[l].size() + g_k_small[rareTarHaps[j]];

      if (g_dash <= f_dash) {
        counter_sel_gf++;
        selectK(htr, common2tot[l - 1], prev_ref_rac_l_com, pbwt_array_A, K,
                tar_hap[htr]);
        length_sel_gf.push(cm_pos[common2tot[l - 1]] - cm_pos[last_reset[htr]]);

        f_dash = 0;
        g_dash = n_ref_haps;
        last_reset[htr] = k - 1;
      }
      ++j;
    } else // common2common
    {
      auto lower0 = std::lower_bound(A_small_idx[l].begin(),
                                     A_small_idx[l].end(), f_k[htr]);
      auto lower1 = std::lower_bound(A_small_idx[l].begin(),
                                     A_small_idx[l].end(), g_k[htr]);
      auto lower2 = std::lower_bound(A_small_idx[l].begin(),
                                     A_small_idx[l].end(), pbwt_index[htr]);

      idx = pbwt_index[htr] - std::distance(A_small_idx[l].begin(), lower2);
      f_dash = f_k[htr] - std::distance(A_small_idx[l].begin(), lower0);
      g_dash = g_k[htr] - std::distance(A_small_idx[l].begin(), lower1);

      if (g_dash <= f_dash) // TODO check this
      {
        counter_sel_gf++;
        if (cm_pos[common2tot[l - 1]] - cm_pos[last_reset[htr]] >
            pbwt_modulo_cm / 5.0f) // only long matches, e.g. 0.02cM
        {
          counter_sel_gf++;
          selectK(htr, common2tot[l - 1], prev_ref_rac_l_com, pbwt_array_A, K,
                  tar_hap[htr]);
          length_sel_gf.push(cm_pos[common2tot[l - 1]] -
                             cm_pos[last_reset[htr]]);
        }
        /*
        selectK(htr, common2tot[l-1], prev_ref_rac_l_com, pbwt_array_A,
        pbwt_depth/4, tar_hap[htr]);
        length_sel_gf.push(cm_pos[common2tot[l-1]]
        - cm_pos[last_reset[htr]]);
        */
        f_dash = 0;
        g_dash = n_ref_haps - A_small_idx[l].size();
        last_reset[htr] = common2tot[l - 1];
      }
    }

    pbwt_index[htr] = idx;
    f_k[htr] = f_dash;
    g_k[htr] = g_dash;
  }
  std::vector<bool> map_big_small(n_ref_haps, true);
  const std::vector<int> &small_idx = A_small_idx[l];
  for (int htr = 0; htr < small_idx.size(); ++htr)
    map_big_small[small_idx[htr]] = false;
  int n_zeros = small_idx[0];
  for (int htr = n_zeros + 1; htr < n_ref_haps; ++htr)
    if (map_big_small[htr])
      pbwt_array_A[n_zeros++] = pbwt_array_A[htr];

  for (int htr = 0; htr < pbwt_small_A.size(); ++htr, ++n_zeros)
    pbwt_array_A[n_zeros] = pbwt_small_A[htr];
}

void haplotype_set::init_rare(const variant_map &M, const int k, const int l) {
  const int A_small_idx_size = A_small_idx[l].size();

  pbwt_small_A.resize(A_small_idx_size);
  for (int i = 0; i < A_small_idx_size; ++i)
    pbwt_small_A[i] = pbwt_array_A[A_small_idx[l][i]];

  pbwt_small_B.resize(A_small_idx_size);
  pbwt_small_V.resize(A_small_idx_size + 1);
  rareTarHaps.clear();

  if (A_small_idx_size <= 0)
    return;

  // Can be precomputed
  for (int kk = k; kk < n_tot_sites; ++kk) {
    if (M.vec_pos[kk]->LQ)
      continue;
    if (flag_common[kk])
      break;
    rareTarHaps.insert(rareTarHaps.end(), SvarTar[kk].begin(),
                       SvarTar[kk].end());
  }

  sort(rareTarHaps.begin(), rareTarHaps.end());
  rareTarHaps.erase(unique(rareTarHaps.begin(), rareTarHaps.end()),
                    rareTarHaps.end());

  for (int id_rare = 0; id_rare < rareTarHaps.size(); ++id_rare) {
    const int htr = rareTarHaps[id_rare];
    int f_dash, g_dash;

    auto lower0 = std::lower_bound(A_small_idx[l].begin(), A_small_idx[l].end(),
                                   f_k[htr]);
    auto lower1 = std::lower_bound(A_small_idx[l].begin(), A_small_idx[l].end(),
                                   g_k[htr]);
    auto lower2 = std::lower_bound(A_small_idx[l].begin(), A_small_idx[l].end(),
                                   pbwt_index[htr]);

    f_dash = std::distance(A_small_idx[l].begin(), lower0);
    g_dash = std::distance(A_small_idx[l].begin(), lower1);

    if (g_dash <= f_dash) {
      // SEL [init rare]
      f_dash = 0;
      g_dash = A_small_idx[l].size();
      last_reset[htr] = k;
    }
    pbwt_small_index[htr] = std::distance(A_small_idx[l].begin(), lower2);
    f_k_small[htr] = f_dash;
    g_k_small[htr] = g_dash;
  }
}

void haplotype_set::selectK(const int htr, const int k, const int ref_rac_l,
                            const std::vector<int> &pbwt_array, const int k0,
                            const unsigned char a) {
  const int sample = tar_hapid2ind[htr];
  const int pbwt_idx = pbwt_index[htr];
  int nh_up = 0;
  int nh_down = 0;
  int d_up = pbwt_index[htr] - f_k[htr];
  int d_down = g_k[htr] - pbwt_index[htr];

  const int k1 = (k0 <= 0) ? K : k0;

  if (d_up < k1 && d_down < k1) {
    nh_up = d_up;
    nh_down = d_down;
  } else if (d_up >= k1 && d_down >= k1) {
    nh_up = k1;
    nh_down = k1;
  } else if (d_up < k1)
    nh_down += k1 - d_up;
  else if (d_down < k1)
    nh_up += k1 - d_down;

  const int max_hh = std::max(nh_up, nh_down);
  for (int o = 0, od = 0, ou = 0; o < max_hh; ++o) {
    if (od < nh_down) {
      const int idx = pbwt_idx + od;
      if (idx < n_ref_haps && ((idx >= ref_rac_l) == a)) {
        pbwt_states[sample][o].push_back(pbwt_array[idx]);
        ++od;
      }
    }
    if (ou < nh_up) {
      const int idx = pbwt_idx - (ou + 1);
      if (idx >= 0 && ((idx >= ref_rac_l) == a)) {
        pbwt_states[sample][o].push_back(pbwt_array[idx]);
        ++ou;
      }
    }
  }
}

void haplotype_set::selectKrare(const int htr, const int k, const int ref_rac_l,
                                const std::vector<int> &pbwt_array,
                                const int k0, const unsigned char a) {
  const int sample = tar_hapid2ind[htr];
  const int pbwt_idx = pbwt_small_index[htr];
  int nh_up = 0;
  int nh_down = 0;
  int d_up = pbwt_small_index[htr] - f_k_small[htr];
  int d_down = g_k_small[htr] - pbwt_small_index[htr];
  const int n_rare_haps = pbwt_array.size();

  const int k1 = (k0 <= 0) ? K : k0;

  if (d_up < k1 && d_down < k1) {
    nh_up = d_up;
    nh_down = d_down;
  } else if (d_up >= k1 && d_down >= k1) {
    nh_up = k1;
    nh_down = k1;
  } else if (d_up < k1)
    nh_down += k1 - d_up;
  else if (d_down < k1)
    nh_up += k1 - d_down;

  const int max_hh = std::max(nh_up, nh_down);
  for (int o = 0, od = 0, ou = 0; o < max_hh; ++o) {
    if (od < nh_down) {
      const int idx = pbwt_idx + od;
      if (idx < n_ref_haps && ((idx >= ref_rac_l) == a))
        pbwt_states[sample][o].push_back(pbwt_array[idx]);
      ++od;
    }
    if (ou < nh_up) {
      const int idx = pbwt_idx - (ou + 1);
      if (idx >= 0 && ((idx >= ref_rac_l) == a))
        pbwt_states[sample][o].push_back(pbwt_array[idx]);
      ++ou;
    }
  }
}

void haplotype_set::performSelection_RARE_INIT_GL(const variant_map &M) {
  if (Kinit == 0) {
    vrb.bullet("No init selection (Kinit=0)");
    return;
  }

  std::vector<int> ref_range_iota(n_ref_haps);
  std::iota(ref_range_iota.begin(), ref_range_iota.end(), 0);
  const int k_init = Kinit;
  const int k_init8 = static_cast<int>(std::floor(k_init * 0.8f));
  std::vector<int> sampled;
  sampled.reserve(k_init);
  std::vector<std::pair<int, int>> tmp_idxHaps_mac;
  tmp_idxHaps_mac.reserve(k_init * 2);

  for (int ind = 0; ind < n_tar_samples; ++ind) {
    const int htr = tar_ind2hapid[ind];
    if (SindTarGL[ind].size() > 0) {
      if (SindTarGL[ind].size() > k_init8) {
        tmp_idxHaps_mac.clear();
        for (auto r0 = 0; r0 < SindTarGL[ind].size(); ++r0)
          tmp_idxHaps_mac.push_back(std::pair<int, int>(
              M.vec_pos[SindTarGL[ind][r0]]->getMAC(), SindTarGL[ind][r0]));
        std::sort(tmp_idxHaps_mac.begin(), tmp_idxHaps_mac.end());

        for (int r0 = 0;
             r0 < SindTarGL[ind].size() && init_states[ind].size() < k_init8;
             r0++) {
          const int idx_rare_variant = tmp_idxHaps_mac[r0].second;
          sampled.clear();
          std::sample(ref_range_iota.begin(),
                      ref_range_iota.begin() + SvarRef[idx_rare_variant].size(),
                      std::back_inserter(sampled), 5, rng.randomEngine);
          for (int i = 0;
               i < sampled.size() && init_states[ind].size() < k_init8; ++i)
            init_states[ind].insert(SvarRef[idx_rare_variant][sampled[i]]);
        }
      } else {
        const int size_max =
            std::max((int)std::floor(k_init8 * 1.0 / SindTarGL[ind].size()), 1);
        for (int r0 = 0; r0 < SindTarGL[ind].size(); r0++) {
          sampled.clear();
          const int idx_rare_variant = SindTarGL[ind][r0];
          std::sample(SvarRef[idx_rare_variant].begin(),
                      SvarRef[idx_rare_variant].end(),
                      std::back_inserter(sampled), size_max, rng.randomEngine);
          for (int r0 = 0;
               r0 < sampled.size() && init_states[ind].size() < k_init8; ++r0)
            init_states[ind].insert(sampled[r0]);
        }
      }
    }
    if (init_states[ind].size() < k_init) {
      sampled.clear();
      std::sample(ref_range_iota.begin(), ref_range_iota.end(),
                  std::back_inserter(sampled), k_init, rng.randomEngine);
      // std::shuffle(std::begin(ref_range_iota), std::end(ref_range_iota),
      // rng.randomEngine);
      for (int r0 = 0; r0 < sampled.size() && init_states[ind].size() < k_init;
           ++r0)
        init_states[ind].insert(sampled[r0]);
    }
  }
  SindTarGL.clear();
  SindTarGL.shrink_to_fit();
}

void haplotype_set::read_list_states(const std::string file_list) {
  input_file fd(file_list);
  if (fd.fail())
    vrb.error("Impossible to open file: " + file_list);

  std::string line;
  int val;
  int hap = 0;
  while (hap < n_tar_haps && getline(fd, line)) {
    std::stringstream s(line);
    while ((s >> val)) {
      if (val < 0 && val >= n_ref_haps)
        vrb.error("State given for haplotype " + std::to_string(hap) +
                  " has an impossible value: " + std::to_string(val));
      list_states[hap].push_back(val);
    }
    if (list_states[hap].size() == 0)
      vrb.error("List of states for haplotype " + std::to_string(hap) +
                " is empty.");
    ++hap;
  }
  if (hap < n_tar_haps)
    vrb.error("Only " + std::to_string(hap) +
              " haplotypes have a list of states. Check the number of rows of "
              "your file.");

  vrb.bullet("List of states read for all haplotypes.");
}
