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
#include <algorithm>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/tmpdir.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <cmath>
#include <containers/haplotype_set.h>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <unordered_set>
#include <objects/genotype.h>
#include <ostream>
#include <random>
#include <set>
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
                                 const int _Kinit, const int _Kpbwt,
                                 const bool _use_mupbwt) {
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
  pbwt_col_bits = std::vector<uint64_t>(n_ref_haps / 64 + 1);
  pbwt_col_prefix = std::vector<int>(n_ref_haps / 64 + 1);
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

  if (Ypacked.size() == 0 && !_use_mupbwt)
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

class mu_score_map {
  static constexpr uint32_t EMPTY = 0xFFFFFFFFu;
  std::vector<uint32_t> keys;
  std::vector<double> vals;
  std::vector<uint32_t> used;
  uint32_t mask;
  uint32_t shift;

  void init(uint32_t cap_log2) {
    keys.assign((size_t)1 << cap_log2, EMPTY);
    vals.assign((size_t)1 << cap_log2, 0.0);
    mask = (1u << cap_log2) - 1;
    shift = 32 - cap_log2;
  }

  void grow() {
    std::vector<uint32_t> old_used;
    old_used.swap(used);
    std::vector<uint32_t> old_keys;
    std::vector<double> old_vals;
    old_keys.swap(keys);
    old_vals.swap(vals);
    init(33 - shift);
    for (uint32_t s : old_used)
      (*this)[old_keys[s]] = old_vals[s];
  }

public:
  mu_score_map() { init(12); }

  double &operator[](uint32_t h) {
    uint32_t s = (h * 2654435761u) >> shift;
    while (keys[s] != EMPTY) {
      if (keys[s] == h)
        return vals[s];
      s = (s + 1) & mask;
    }
    if ((used.size() + 1) * 2 > keys.size()) {
      grow();
      return (*this)[h];
    }
    keys[s] = h;
    used.push_back(s);
    return vals[s];
  }

  void clear() {
    for (uint32_t s : used) {
      keys[s] = EMPTY;
      vals[s] = 0.0;
    }
    used.clear();
  }
};

static size_t mu_batch_size(size_t n, int threads, int cap) {
  size_t b = n / ((size_t)std::max(threads, 1) * 4);
  return std::max<size_t>(1, std::min<size_t>(b, (size_t)std::max(cap, 1)));
}

void haplotype_set::matchHapsFromMuPBWT(
    pbwt &mupbwt, const variant_map &V, const bool main_iteration,
    const genotype_set &G, uint32_t n_sites, int threads, haplotype_set &H,
    float min_cm, float short_cm, float medium_cm, int ma, int mc, int md,
    bool persistence_enabled, double persistence_decay,
    double persistence_floor, int batch_cap) {
  omp_set_num_threads(threads);

  if (Kpbwt == 0 || Kpbwt >= n_ref_haps) {
    vrb.bullet("No PBWT selection (Kpbwt=0 or Kpbwt >= n_ref_haps)");
    return;
  }

  assert(cm_pos.size() == n_sites);

  tac.clock();

  for (int e = 0; e < n_tar_samples; e++) {
    for (int j = 0; j < K; ++j) {
      std::vector<int>().swap(pbwt_states[e][j]);
    }
  }

  size_t chunk_size = mc;
  auto tmp_d = pbwt_depth;
  pbwt_depth = md;

  const double thr_min_cm = min_cm;
  const double thr_short_cm = short_cm;
  const double thr_medium_cm = medium_cm;
  const double region_cm = std::max(1e-6, (double)(cm_pos.back() - cm_pos.front()));

  uint32_t n_haps = mupbwt.n_haps;

  uint32_t n = G.vecG.size() * 2;
  if (prev_selected.size() != n)
    prev_selected.assign(n, {});
  const uint32_t MAX_EXACT_STEPS = K * ma;
  double sum_smem_time = 0.0;
  double sum_rank_time = 0.0;
  std::vector<uint8_t> site_minor;
  if (mu_scan_view.off.empty()) {
    pbwt_build_scan_view(&mupbwt, mu_scan_view, true);
  }
  site_minor.assign(n_sites, 2);
  for (size_t s = 0; s < n_sites; ++s)
    if (!flag_common[s])
      site_minor[s] = major_alleles[s] ? 0 : 1;
  size_t diag_sel_sum = 0;
  size_t diag_sel_max = 0;
  size_t diag_pool_sum = 0;
  size_t diag_pool_max = 0;
  struct mu_query_ctx {
    mu_score_map haplo_score;
    std::vector<uint32_t> active_haplos;
    std::vector<std::vector<int>> local_states;
    std::vector<bool> rescued_ids;
    std::vector<int> local_rescued;
    std::vector<uint8_t> query;
    std::vector<int> q_rare;
    size_t q;
    uint32_t p, p_p, l, p_l;
    uint32_t c_p, c_i, c_r, p_i;
    uint8_t c_s, p_s;
  };
  const size_t n_batch = mu_batch_size(n, threads, batch_cap);
  const size_t n_blocks = (n + n_batch - 1) / n_batch;
#pragma omp parallel num_threads(threads) default(none)                        \
    shared(mupbwt, G, n, n_sites, chunk_size, pbwt_depth, pbwt_states,         \
               tar_hapid2ind, thr_min_cm, thr_short_cm, thr_medium_cm,         \
               region_cm, n_haps, K,       \
               sum_smem_time, sum_rank_time, MAX_EXACT_STEPS, H,               \
               diag_sel_sum, diag_sel_max, diag_pool_sum,     \
               diag_pool_max, prev_selected, persistence_enabled,              \
               persistence_decay, persistence_floor, site_minor, n_batch,      \
               n_blocks)
  {
    std::vector<mu_query_ctx> ctx(std::min(n_batch, (size_t)n));
    for (auto &Q : ctx) {
      Q.active_haplos.reserve(2000);
      Q.local_states.resize(pbwt_depth);
      Q.rescued_ids.assign(n_haps, false);
      Q.query.assign(n_sites, 0);
    }

    pbwt_col *all_cols = mupbwt.cols.a;
    const bool use_view = !mu_scan_view.off.empty();
    const uint32_t *v_off = mu_scan_view.off.data();
    const uint32_t *v_pool = mu_scan_view.pool.data();

    auto get_l_ref = [&](mu_query_ctx &Q, uint32_t end_site,
                         uint32_t h) -> uint32_t {
      const std::vector<int> &q_rare = Q.q_rare;
      const uint8_t *query = Q.query.data();
      const int e = (int)end_site;
      const std::vector<int> &sh = ShapRef[h];
      int a = (int)(std::upper_bound(q_rare.begin(), q_rare.end(), e) - q_rare.begin()) - 1;
      int b = (int)(std::upper_bound(sh.begin(), sh.end(), e) - sh.begin()) - 1;
      while (a >= 0 && b >= 0 && q_rare[a] == sh[b]) {
        a--;
        b--;
      }
      int mis = std::max(a >= 0 ? q_rare[a] : -1, b >= 0 ? sh[b] : -1);
      int c = (int)(std::upper_bound(common2tot.begin(), common2tot.end(), e) - common2tot.begin()) - 1;
      const size_t row_bytes = HvarRef.n_cols >> 3;
      const unsigned char *hcol = HvarRef.bytes + (h >> 3);
      for (; c >= 0 && common2tot[c] > mis; c--) {
        if (c >= 8)
          __builtin_prefetch(hcol + (size_t)(c - 8) * row_bytes, 0, 0);
        if (HvarRef.get(c, h) != query[common2tot[c]]) {
          mis = common2tot[c];
          break;
        }
      }
      return (uint32_t)(e - mis);
    };
    auto compute_match_cm = [&](uint32_t end_col, uint32_t len_sites) -> double {
      if (len_sites == 0)
        return 0.0;
      uint32_t start_col = (end_col + 1 > len_sites) ? (end_col + 1 - len_sites) : 0;
      return std::max(0.0f, cm_pos[end_col] - cm_pos[start_col]);
    };
    auto apply_score_inline = [&](mu_query_ctx &Q, uint32_t start_p,
                                  uint32_t match_len, uint32_t col,
                                  double match_cm, bool last = false) {
      mu_score_map &haplo_score = Q.haplo_score;
      std::vector<uint32_t> &active_haplos = Q.active_haplos;
      double contrib = match_cm / region_cm;
      double weight = contrib;
      uint32_t depth_limit = 0;
      bool exact_only = false;

      if (last) {
        exact_only = true;
        depth_limit = K * 2;
      } else {
        if (match_cm < thr_short_cm) {
          depth_limit = K;
          exact_only = true;
        } else if (match_cm < thr_medium_cm) {
          depth_limit = K * 2;
          exact_only = true;
        } else {
          exact_only = true;
          depth_limit = MAX_EXACT_STEPS;
        }
      }

      auto add_score = [&](uint32_t h_id) {
        if (haplo_score[h_id] == 0.0) {
          active_haplos.push_back(h_id);
        }
        haplo_score[h_id] += weight;
      };

      add_score(start_p);
      if (last && !Q.rescued_ids[start_p]) {
        Q.rescued_ids[start_p] = true;
        Q.local_rescued.push_back(start_p);
      }
      uint32_t t_p = start_p;
      uint32_t steps = 0;
      while (steps < depth_limit) {
        uint32_t d_p = phi_inv_f(&mupbwt.phid, t_p, col + 1);
        if (d_p == n_haps)
          break;
        if (exact_only) {
          uint32_t t_l = phi_l_fast(&mupbwt.phid, d_p, col + 1);
          if (t_l >= match_len) {
            add_score(d_p);
            t_p = d_p;
          } else if (last && t_l >= match_len / 4) {
            add_score(d_p);
            t_p = d_p;
          } else {
            break;
          }
        } else {
          add_score(d_p);
          t_p = d_p;
        }
        steps++;
      }

      t_p = start_p;
      steps = 0;
      while (steps < depth_limit) {
        uint32_t t_l;
        uint32_t u_p = phi_f_l(&mupbwt.phid, t_p, col + 1, &t_l);
        if (u_p == n_haps)
          break;
        if (exact_only) {
          if (t_l >= match_len) {
            add_score(u_p);
            t_p = u_p;
          } else if (last && t_l >= match_len / 4) {
            add_score(u_p);
            t_p = u_p;
          } else {
            break;
          }
        } else {
          add_score(u_p);
          t_p = u_p;
        }
        steps++;
      }
    };

    auto scan_site = [&](mu_query_ctx &Q, size_t i, pbwt_col *c_col,
                         pbwt_col *next_col) {
      const uint8_t *query = Q.query.data();
      uint8_t q_s = query[i];
      if (LIKELY(q_s == Q.c_s)) {
        Q.p = Q.c_p;
        Q.l = (i != 0 && Q.p_p == Q.p) ? Q.p_l + 1 : 1;
      } else {
        uint32_t t = c_arr_get(&c_col->t, Q.c_r);
        uint32_t p_n = c_col->p.n;

        if (UNLIKELY(p_n == 1)) {
          Q.p = n_haps;
          Q.l = 0;
          Q.p_s = 0;
          if (next_col) {
            Q.c_p = c_arr_get(&next_col->e_pa, next_col->e_pa.n - 1);
            Q.c_i = n_haps - 1;
            Q.c_r = next_col->e_pa.n - 1;
            Q.c_s = get_ns(next_col->zero, Q.c_r);
          }
          goto check_match;
        }

        uint32_t d;
        bool is_up = (Q.c_r != 0 && ((Q.c_i < t) || (Q.c_r == p_n - 1)));

        if (is_up) {
          Q.p_i = Q.c_i;
          Q.c_i = c_arr_get(&c_col->p, Q.c_r) - 1;
          Q.c_p = c_arr_get(&c_col->e_pa, Q.c_r - 1);
          Q.c_r--;
          d = Q.p_i - Q.c_i;
        } else {
          Q.p_i = Q.c_i;
          Q.c_i = c_arr_get(&c_col->p, Q.c_r + 1);
          Q.c_p = c_arr_get(&c_col->b_pa, Q.c_r + 1);
          Q.c_r++;
          d = Q.c_i - Q.p_i;
        }

        if (d < 30) {
          uint32_t m = is_up ? get_l_u(&mupbwt, Q.p_p, d, i)
                             : get_l_d(&mupbwt, Q.p_p, d, i);
          Q.l = std::min((uint32_t)m, Q.p_l) + 1;
        } else {
          Q.l = get_l_ref(Q, i, Q.c_p);
        }
        Q.p = Q.c_p;
      }
      Q.p_s = Q.c_s;
      if (next_col) {
        if (use_view) {
          const uint32_t *cv = v_pool + v_off[i];
          Q.c_i = scan_hdr_n(cv[0]) ? Q.c_i + cv[1 + Q.c_r] : fl(c_col, Q.c_i, Q.c_r);
          const uint32_t *nv = v_pool + v_off[i + 1];
          const uint32_t nv_n = scan_hdr_n(nv[0]);
          if (nv_n) {
            const uint32_t *st = nv + 1 + nv_n;
            uint32_t r = 0;
            for (uint32_t k = 0; k + 1 < nv_n; k++)
              r += Q.c_i >= st[k];
            Q.c_r = r;
          } else if (scan_hdr_samples(nv[0])) {
            Q.c_r = get_r_sampled(next_col, nv, Q.c_i);
          } else {
            Q.c_r = get_r(next_col, Q.c_i);
          }
          Q.c_s = scan_hdr_ns0(nv[0]) ^ (Q.c_r & 1);
        } else {
          Q.c_i = fl(c_col, Q.c_i, Q.c_r);
          Q.c_r = get_r(next_col, Q.c_i);
          Q.c_s = get_ns(next_col->zero, Q.c_r);
        }
      }

    check_match:
      if (i > 0 && Q.p_p != Q.p && Q.p_l > 0 && Q.p_l >= Q.l) {
        double match_cm_1 = compute_match_cm(i - 1, Q.p_l);
        if (match_cm_1 > thr_min_cm)
          apply_score_inline(Q, Q.p_p, Q.p_l, i - 1, match_cm_1);
      }

      if (i == n_sites - 1 && Q.l != 0) {
        double match_cm_2 = compute_match_cm(i, Q.l);
        if (match_cm_2 > thr_min_cm)
          apply_score_inline(Q, Q.p, Q.l, i, match_cm_2);
      }

      if (Q.p_s == 1 && (use_view ? scan_hdr_rare15(v_pool[v_off[i]])
                                  : c_col->e_pa.n <= 15)) {
        double match_cm_3 = compute_match_cm(i, Q.l);
        apply_score_inline(Q, Q.p, Q.l, i, match_cm_3, true);
      }

      Q.p_p = Q.p;
      Q.p_l = Q.l;
    };

#pragma omp for schedule(dynamic)
    for (size_t blk = 0; blk < n_blocks; ++blk) {
      const size_t q_begin = blk * n_batch;
      const size_t n_q = std::min(n_batch, (size_t)n - q_begin);

      for (size_t b = 0; b < n_q; ++b) {
        mu_query_ctx &Q = ctx[b];
        Q.q = q_begin + b;
        Q.haplo_score.clear();
        for (uint32_t h : Q.active_haplos)
          Q.rescued_ids[h] = false;
        Q.active_haplos.clear();
        Q.local_rescued.clear();
        for (auto &vec : Q.local_states) {
          vec.clear();
        }

        size_t ind_idx = Q.q / 2;
        bool is_h1 = (Q.q % 2 != 0);
        const auto &genotype_ptr = G.vecG[ind_idx];
        const auto &source_H = is_h1 ? genotype_ptr->H1 : genotype_ptr->H0;
        uint8_t *query = Q.query.data();

        for (size_t s = 0; s < n_sites; ++s) {
          query[s] = source_H[s] ? 1 : 0;
        }

        Q.q_rare.clear();
        for (size_t s = 0; s < n_sites; ++s)
          if (site_minor[s] == query[s])
            Q.q_rare.push_back((int)s);
        Q.p = 0;
        Q.p_p = 0;
        Q.l = 0;
        Q.p_l = 0;
        pbwt_col *c_col = &all_cols[0];
        Q.c_p = c_arr_get(&c_col->e_pa, c_col->e_pa.n - 1);
        Q.c_i = n_haps - 1;
        Q.c_r = c_col->e_pa.n - 1;
        Q.c_s = get_ns(c_col->zero, Q.c_r);
        Q.p_i = Q.c_i;
        Q.p_s = Q.c_s;
      }

      for (size_t i = 0; i < n_sites; i++) {
        pbwt_col *c_col = &all_cols[i];
        pbwt_col *next_col = (i < n_sites - 1) ? &all_cols[i + 1] : NULL;
        if (use_view) {
          if (i + PBWT_PF_STRUCT < n_sites &&
              scan_hdr_n(v_pool[v_off[i + PBWT_PF_STRUCT]]) == 0) {
            const char *a = (const char *)&all_cols[i + PBWT_PF_STRUCT];
            __builtin_prefetch(a, 0, 3);
            __builtin_prefetch(a + 64, 0, 3);
            __builtin_prefetch(a + sizeof(pbwt_col) - 1, 0, 3);
          }
          if (i + PBWT_PF_DATA < n_sites &&
              scan_hdr_n(v_pool[v_off[i + PBWT_PF_DATA]]) == 0) {
            const uint32_t *pv = v_pool + v_off[i + PBWT_PF_DATA];
            const uint32_t pns = scan_hdr_samples(pv[0]);
            if (pns > 0 && pns <= 128)
              for (uint32_t k = 0; k <= pns; k += 16)
                __builtin_prefetch(pv + k, 0, 3);
            __builtin_prefetch(all_cols[i + PBWT_PF_DATA].p.data, 0, 3);
            __builtin_prefetch(all_cols[i + PBWT_PF_DATA].uv.data, 0, 3);
          }
        } else {
          pbwt_prefetch_cols(all_cols, i, n_sites);
        }
        for (size_t b = 0; b < n_q; ++b)
          scan_site(ctx[b], i, c_col, next_col);
      }

      for (size_t b = 0; b < n_q; ++b) {
        mu_query_ctx &Q = ctx[b];
        const size_t q = Q.q;
        uint32_t target_ind = tar_hapid2ind[q];
        mu_score_map &haplo_score = Q.haplo_score;
        std::vector<uint32_t> &active_haplos = Q.active_haplos;
        std::vector<std::vector<int>> &local_states = Q.local_states;

        if (persistence_enabled) {
          for (const auto &pr : prev_selected[q]) {
            int h = pr.first;
            double decayed = pr.second * persistence_decay;
            if (decayed < persistence_floor)
              continue;
            if (haplo_score[h] == 0.0) {
              active_haplos.push_back((uint32_t)h);
              haplo_score[h] = decayed;
            }
          }
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
        size_t elements_needed = (pbwt_depth)*chunk_size;
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
          uint32_t t_p = Q.p;
          uint32_t steps = 0;
          local_states[pbwt_depth - 1].push_back(Q.p);

          while (steps < chunk_size) {
            uint32_t d_p = phi_inv_f(&mupbwt.phid, t_p, n_sites);
            if (d_p == n_haps)
              break;
            local_states[pbwt_depth - 1].push_back(d_p);
            t_p = d_p;
            steps++;
          }
        }

        size_t diag_local_total = 0;
        for (size_t d = 0; d < pbwt_depth; d++) {
          diag_local_total += local_states[d].size();
        }

        if (persistence_enabled) {
          std::vector<std::pair<int, double>> next_prev;
          next_prev.reserve(diag_local_total);
          for (size_t d = 0; d < (size_t)pbwt_depth; d++)
            for (int h : local_states[d])
              next_prev.emplace_back(h, std::min(haplo_score[h], 1.0));
          prev_selected[q] = std::move(next_prev);
        }

#pragma omp critical
        {
          pbwt_states[target_ind][0].insert(pbwt_states[target_ind][0].end(),
                                            Q.local_rescued.begin(),
                                            Q.local_rescued.end());
          for (size_t d = 0; d < pbwt_depth; d++) {
            if (!local_states[d].empty()) {
              pbwt_states[target_ind][d].insert(pbwt_states[target_ind][d].end(),
                                                local_states[d].begin(),
                                                local_states[d].end());
            }
          }
          diag_sel_sum += diag_local_total;
          diag_sel_max = std::max(diag_sel_max, diag_local_total);
          diag_pool_sum += haplo_vec.size();
          diag_pool_max = std::max(diag_pool_max, haplo_vec.size());
        }
      }
    }
  }

  vrb.bullet("Mu-PBWT selection (" + stb.str(tac.rel_time() * 1.0 / 1000, 2) +
             "s)");
  pbwt_depth = tmp_d;
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
  uint64_t *bits = pbwt_col_bits.data();
  int *prefix = pbwt_col_prefix.data();
  std::fill(pbwt_col_bits.begin(), pbwt_col_bits.end(), 0);
  prefix[0] = 0;
  while (m < size_a) {
    z = *pY++;              // read the encoded symbol
    n = p3decode[z & 0x7f]; // get how many chars are encoded
    mm = m + n;             // get limit of this block
    z >>= 7;                // get the symbol itself
    std::copy(pbwt_array_A.begin() + m, pbwt_array_A.begin() + mm,
              pbwt_array_B.begin() + *occ[z]); // update B, (next A)
    for (int w = (m >> 6) + 1; (w << 6) <= mm; ++w)
      prefix[w] = v - ref_rac_l + z * ((w << 6) - m);
    if (z) {
      const int w0 = m >> 6, w1 = (mm - 1) >> 6;
      const uint64_t mask0 = ~0ULL << (m & 63),
                     mask1 = ~0ULL >> (63 - ((mm - 1) & 63));
      if (w0 == w1)
        bits[w0] |= mask0 & mask1;
      else {
        bits[w0] |= mask0;
        std::fill(bits + w0 + 1, bits + w1, ~0ULL);
        bits[w1] |= mask1;
      }
    }
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

    const int v_idx = pbwt_rank(pbwt_index[htr]);
    const int v_f = pbwt_rank(f_k[htr]);
    const int v_g = pbwt_rank(g_k[htr]);
    idx = tar_hap[htr] * (ref_rac_l + v_idx) +
          (1 - tar_hap[htr]) * (pbwt_index[htr] - v_idx);
    f_dash = tar_hap[htr] * (ref_rac_l + v_f) +
             (1 - tar_hap[htr]) * (f_k[htr] - v_f);
    g_dash = tar_hap[htr] * (ref_rac_l + v_g) +
             (1 - tar_hap[htr]) * (g_k[htr] - v_g);

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
  const std::vector<int> &small_idx = A_small_idx[l];
  const int n_small = small_idx.size();
  int *A = pbwt_array_A.data();
  int n_kept = small_idx[0];
  for (int i = 0; i < n_small; ++i) {
    const int gap_start = small_idx[i] + 1;
    const int gap_end = (i + 1 < n_small) ? small_idx[i + 1] : n_ref_haps;
    std::copy(A + gap_start, A + gap_end, A + n_kept);
    n_kept += gap_end - gap_start;
  }
  std::copy(pbwt_small_A.begin(), pbwt_small_A.end(), A + n_kept);
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
