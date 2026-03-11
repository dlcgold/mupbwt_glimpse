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

#include <caller/caller_header.h>
#include <omp.h>
void *phase_callback(void *ptr) {
  caller *S = static_cast<caller *>(ptr);
  int id_worker, id_job;
  float prog_step = 1.0 / S->G.n_ind;
  const int n_ind = S->G.n_ind;
  pthread_mutex_lock(&S->mutex_workers);
  id_worker = S->i_workers++;
  pthread_mutex_unlock(&S->mutex_workers);
  for (;;) {
    pthread_mutex_lock(&S->mutex_workers);
    id_job = S->i_jobs++;
    if (id_job <= n_ind)
      vrb.progress("  * HMM imputation", id_job * prog_step);
    pthread_mutex_unlock(&S->mutex_workers);
    if (id_job < n_ind)
      S->phase_individual(id_worker, id_job);
    else
      pthread_exit(NULL);
  }
  return nullptr;
}

void caller::phase_individual(const int id_worker, const int id_job) {
  const int ploidy = G.vecG[id_job]->ploidy;

  COND[id_worker]->select(id_job, current_stage);

  if (current_stage == STAGE_INIT)
    G.vecG[id_job]->initHaplotypeLikelihoods(HLC[id_worker], min_gl);
  else {
    if (ploidy > 1)
      G.vecG[id_job]->makeHaplotypeLikelihoods(HLC[id_worker], true, min_gl);
    else
      G.vecG[id_job]->initHaplotypeLikelihoods(HLC[id_worker], min_gl);
  }
  HMM[id_worker]->computePosteriors(HLC[id_worker], G.vecG[id_job]->flat,
                                    HP0[id_worker]);
  G.vecG[id_job]->sampleHaplotypeH0(HP0[id_worker]);
  if (ploidy > 1) {
    G.vecG[id_job]->makeHaplotypeLikelihoods(HLC[id_worker], false, min_gl);
    HMM[id_worker]->computePosteriors(HLC[id_worker], G.vecG[id_job]->flat,
                                      HP1[id_worker]);
    G.vecG[id_job]->sampleHaplotypeH1(HP1[id_worker]);
    DMM[id_worker]->rephaseHaplotypes(G.vecG[id_job]->H0, G.vecG[id_job]->H1,
                                      G.vecG[id_job]->flat);
  }

  if (current_stage == STAGE_MAIN) {
    if (ploidy > 1)
      G.vecG[id_job]->storeGenotypePosteriorsAndHaplotypes(HP0[id_worker],
                                                           HP1[id_worker]);
    else
      G.vecG[id_job]->storeGenotypePosteriorsAndHaplotypes(HP0[id_worker]);
  }

  if (options["threads"].as<int>() > 1)
    pthread_mutex_lock(&mutex_workers);
  statH.push(COND[id_worker]->n_states * 1.0);
  statC.push(COND[id_worker]->polymorphic_sites.size() * 100.0 /
             COND[id_worker]->n_tot_sites);
  if (options["threads"].as<int>() > 1)
    pthread_mutex_unlock(&mutex_workers);
}

void caller::phase_iteration() {
  tac.clock();
  int n_thread = options["threads"].as<int>();
  i_workers = 0;
  i_jobs = 0;
  statH.clear();
  statC.clear();

  float prog_step = 1.0 / G.n_ind;
  float prog_bar = 0.0;

  if (n_thread > 1) {
    for (int t = 0; t < n_thread; t++)
      pthread_create(&id_workers[t], NULL, phase_callback,
                     static_cast<void *>(this));
    for (int t = 0; t < n_thread; t++)
      pthread_join(id_workers[t], NULL);
  } else
    for (int i = 0; i < G.n_ind; i++) {
      phase_individual(0, i);
      prog_bar += prog_step;
      vrb.progress("  * HMM imputation", prog_bar);
    }
  vrb.bullet("HMM imputation [#states=" + stb.str(statH.mean(), 1) +
             " / %poly=" + stb.str(statC.mean(), 1) + "%] (" +
             stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");
}

void caller::phase_loop() {
  if (use_mu) {
    const char *env = std::getenv("OMP_NUM_THREADS");
    if (env) {
      omp_set_num_threads(std::atoi(env));
    } else {
      omp_set_num_threads(options["threads"].as<int>());
    }
  }
  // First Iteration
  current_stage = STAGE_INIT;
  vrb.title("Initializing iteration");

  H.initRareTar(G, V);
  H.performSelection_RARE_INIT_GL(V); // not parallel

  phase_iteration();

  H.init_states.clear();
  H.init_states.shrink_to_fit();

  // Burn-in 0
  current_stage = STAGE_BURN;
  int nBurnin = options["burnin"].as<int>();
  // nBurnin = 1;
  //    std::vector<std::string> tmp;
  //    rlpbwt_int mu(options["reference_panel"].as<std::string>(), tmp, H, V,
  //                  V.input_gregion, 1);
  // mupbwt = mu;
  std::vector<int> sites;

  vrb.bullet("mu-PBWT: " + stb.str(mupbwt.n_haps) + " haplotypes and " +
             stb.str(mupbwt.n_sites) + " variants");
  for (auto a : H.pbwt_grp) {
    sites.push_back(H.common2tot[a]);
    // std::cout << H.common2tot[a] << "\t";
  }
  //    for (int i = 0; i < mupbwt.height; i++) {
  //        std::cerr << mupbwt.get_row(i) << "\n";
  //    }

  int n_q = 0;
  vrb.bullet("PBWT sites " + stb.str(sites.size()));

  vrb.bullet("# target " + stb.str(G.vecG.size() * 2));
  for (int iter = 0; iter < nBurnin; iter++) {
    vrb.title("Burn-in iteration [" + stb.str(iter + 1) + "/" +
              stb.str(nBurnin) + "]");

    H.updateHaplotypes(G);
    H.transposeRareTar();
    if (use_mu) {
      tac.clock();
      uint32_t n = G.vecG.size() * 2;
      uint32_t n_sites = mupbwt.n_sites;

      uint8_t *q_p = nullptr;
      size_t mat_size = (size_t)n * n_sites * sizeof(uint8_t);

      if (posix_memalign((void **)&q_p, 64, mat_size) != 0) {
        fprintf(stderr, "Memory allocation error for q_p\n");
        exit(-1);
      }
      memset(q_p, 0, mat_size);

      uint32_t q_idx = 0;

      for (auto &i : G.vecG) {
        uint32_t site_c = 0;
        uint32_t pbwt_site = 0;

        for (auto &&j : i->H0) {
          if (!use_mu_common || H.flag_common[site_c]) {
            if (j) {
              q_p[(size_t)q_idx * n_sites + pbwt_site] = 1;
            }
            pbwt_site++;
          }
          site_c++;
        }
        q_idx++;

        site_c = 0;
        pbwt_site = 0;

        for (auto &&j : i->H1) {
          if (!use_mu_common || H.flag_common[site_c]) {
            if (j) {
              q_p[(size_t)q_idx * n_sites + pbwt_site] = 1;
            }
            pbwt_site++;
          }
          site_c++;
        }
        q_idx++;
      }
      vrb.bullet("Mu-PBWT loaded queries (" +
                 stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");
      // if (!use_smems && !use_mpsc) {
      //   H.matchHapsFromMuPBWT(mupbwt, V, false, sites, queries);
      // } else if (use_smems && !use_mpsc) {
      //   H.matchHapsFromMuPBWTSMEMS(mupbwt, V, false, sites, queries);
      // } else if (!use_smems && use_mpsc) {
      //   H.matchHapsFromMuPBWTSMEMS(mupbwt, V, false, sites, queries);
      // }
      H.matchHapsFromMuPBWT(mupbwt, V, false, sites, q_p, n, n_sites,
                            options["threads"].as<int>());
      n_q = n;
      free(q_p);
    } else {
      H.matchHapsFromCompressedPBWTSmall(V, false);
    }
    // for (int i = 0; i < H.pbwt_states[0].size(); i++) {
    //   vrb.bullet("after extraction pbwt_states: " + stb.str(n_q) +
    //              " queries, " + stb.str(i) + " depth, " +
    //              stb.str(H.pbwt_states[0][i].size()) + " haplotypes");
    // }
    // vrb.bullet("n_q = " + std::to_string(n_q));
    // vrb.bullet("pbwt_states.size() = " +
    // std::to_string(H.pbwt_states.size()));
    size_t max_depth = 0;
    for (size_t q = 0; q < n_q / 2; q++) {
      max_depth = std::max(max_depth, H.pbwt_states[q].size());
    }

    for (size_t d = 0; d < max_depth; d++) {
      size_t total_haplos = 0;

      for (size_t q = 0; q < n_q / 2; q++) {
        if (H.pbwt_states[q].size() > d) {
          total_haplos += H.pbwt_states[q][d].size();
        }
      }

      vrb.bullet("after extraction pbwt_states: total haplotypes at depth " +
                 std::to_string(d) + " = " + std::to_string(total_haplos));
    }
    current_stage = STAGE_RESTRICT;
    phase_iteration();
    current_stage = STAGE_BURN;
    phase_iteration();
  }
  vrb.bullet("having " + stb.str(sites.size()) + " sites");
  // Main
  current_stage = STAGE_MAIN;
  int nMain = options["main"].as<int>();
  // nMain = 0;
  for (int iter = 0; iter < nMain; iter++) {
    vrb.title("Main iteration [" + stb.str(iter + 1) + "/" + stb.str(nMain) +
              "]");
    H.updateHaplotypes(G);
    H.transposeRareTar();

    if (use_mu) {
      tac.clock();
      uint32_t n = G.vecG.size() * 2;
      uint32_t n_sites = mupbwt.n_sites;

      uint8_t *q_p = nullptr;
      size_t mat_size = (size_t)n * n_sites * sizeof(uint8_t);

      if (posix_memalign((void **)&q_p, 64, mat_size) != 0) {
        fprintf(stderr, "Memory allocation error for q_p\n");
        exit(-1);
      }
      memset(q_p, 0, mat_size);

      uint32_t q_idx = 0;

      for (auto &i : G.vecG) {
        uint32_t site_c = 0;
        uint32_t pbwt_site = 0;

        for (auto &&j : i->H0) {
          if (!use_mu_common || H.flag_common[site_c]) {
            if (j) {
              q_p[(size_t)q_idx * n_sites + pbwt_site] = 1;
            }
            pbwt_site++;
          }
          site_c++;
        }
        q_idx++;

        site_c = 0;
        pbwt_site = 0;

        for (auto &&j : i->H1) {
          if (!use_mu_common || H.flag_common[site_c]) {
            if (j) {
              q_p[(size_t)q_idx * n_sites + pbwt_site] = 1;
            }
            pbwt_site++;
          }
          site_c++;
        }
        q_idx++;
      }
      vrb.bullet("Mu-PBWT loaded queries (" +
                 stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");
      H.matchHapsFromMuPBWT(mupbwt, V, true, sites, q_p, n, n_sites,
                            options["threads"].as<int>());
      // if (!use_smems && !use_mpsc) {
      //   H.matchHapsFromMuPBWT(mupbwt, V, false, sites, queries);
      // } else if (use_smems && !use_mpsc) {
      //   H.matchHapsFromMuPBWTSMEMS(mupbwt, V, false, sites, queries);
      // } else if (!use_smems && use_mpsc) {
      //   H.matchHapsFromMuPBWTMPSC(mupbwt, V, false, sites, queries);
      // }
      n_q = n;
      free(q_p);
    } else {
      H.matchHapsFromCompressedPBWTSmall(V, true);
    }
    size_t max_depth = 0;
    for (size_t q = 0; q < n_q / 2; q++) {
      max_depth = std::max(max_depth, H.pbwt_states[q].size());
    }

    for (size_t d = 0; d < max_depth; d++) {
      size_t total_haplos = 0;

      for (size_t q = 0; q < n_q / 2; q++) {
        if (H.pbwt_states[q].size() > d) {
          total_haplos += H.pbwt_states[q][d].size();
        }
      }

      vrb.bullet("after extraction pbwt_states: total haplotypes at depth " +
                 std::to_string(d) + " = " + std::to_string(total_haplos));
    }
    phase_iteration();
  }

  // Finalization
  // vrb.title("Finalization");
  for (int i = 0; i < G.vecG.size(); i++)
    G.vecG[i]->sortAndNormAndInferGenotype();

  // vrb.bullet("done");
}
