/* Copyright (C) 2015-16 University of Southern California and
 *                       Andrew D. Smith and Jenny Qu
 *
 * Authors: Jenny Qu and Andrew Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "epiphy_mcmc.hpp"
#include "epiphy_utils.hpp"
#include "sufficient_statistics_helpers.hpp"

#include <vector>
#include <random>
#include <iostream>
#include <cmath>       /* pow */
#include <algorithm>  //std::max, min
#include <boost/math/distributions/students_t.hpp>

using std::vector;
using std::cerr;
using std::endl;
using std::max;


#include <functional>
using std::placeholders::_1;
using std::bind;
using std::plus;

static const double PROBABILITY_GUARD = 1e-10;
static const double CBM_THETA = 0.5;
static const double CBM_EPS = 1e-4;


void
sum(const vector<mcmc_stat> &mcmcstats,
    mcmc_stat &sum_mcmc_stat) {
  sum_mcmc_stat = mcmcstats[0];
  for (size_t i = 1; i < mcmcstats.size(); ++i) {
    sum_mcmc_stat.root_start_distr.first += mcmcstats[i].root_start_distr.first;
    sum_mcmc_stat.root_start_distr.second += mcmcstats[i].root_start_distr.second;
    sum_mcmc_stat.root_distr += mcmcstats[i].root_distr;
    for (size_t j = 0; j < mcmcstats[0].start_distr.size(); ++j) {
      sum_mcmc_stat.start_distr[j] += mcmcstats[i].start_distr[j];
      sum_mcmc_stat.triad_distr[j] += mcmcstats[i].triad_distr[j];
    }
  }
}

void
average(const vector<mcmc_stat> &mcmcstats,
        mcmc_stat &ave_mcmc_stat) {
  sum(mcmcstats, ave_mcmc_stat);
  const size_t N = mcmcstats.size();
  ave_mcmc_stat.root_start_distr.first /= N;
  ave_mcmc_stat.root_start_distr.second /= N;
  ave_mcmc_stat.root_distr.div(N);
  for (size_t j = 0; j < mcmcstats[0].start_distr.size(); ++j) {
    ave_mcmc_stat.start_distr[j].div(N);
    ave_mcmc_stat.triad_distr[j].div(N);
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////   MCMC output analysis      //////////////////////////////////
/// multiple chain , EPSR (estimated potential scale reduction factor) /////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void
mcmc_stat::scale() {
  double tot = root_start_distr.first + root_start_distr.second;
  assert(tot > 0);
  root_start_distr.first = root_start_distr.first/tot;
  root_start_distr.second = root_start_distr.second/tot;

  tot = root_distr.uu + root_distr.um + root_distr.mu + root_distr.mm;
  assert(tot > 0);
  root_distr.div(tot);

  assert(start_distr.size() > 1 && triad_distr.size() > 1);
  for (size_t i = 1; i < start_distr.size(); ++i) {
    start_distr[i].to_probabilities();
    tot = (triad_distr[i].uuu + triad_distr[i].uum +
           triad_distr[i].umu + triad_distr[i].umm +
           triad_distr[i].muu + triad_distr[i].mum +
           triad_distr[i].mmu + triad_distr[i].mmm);
    assert(tot > 0);
    triad_distr[i].div(tot);
  }
}


static void
flatten_mcmcstat(const mcmc_stat &m,
                 vector<double> &v) {
  v.clear();
  v.push_back(m.root_start_distr.first);
  v.push_back(m.root_start_distr.second);

  vector<double> p;
  m.root_distr.flatten(p);
  v.insert(v.end(), p.begin(), p.end());

  assert(m.start_distr.size() > 1 && m.triad_distr.size() > 1);
  for (size_t i = 1; i < m.start_distr.size(); ++i) {
    m.start_distr[i].flatten(p);
    v.insert(v.end(), p.begin(), p.end());
  }

  for (size_t i = 1; i < m.triad_distr.size(); ++i) {
    m.triad_distr[i].flatten(p);
    v.insert(v.end(), p.begin(), p.end());
  }
}


void
within_mean(const vector<vector<vector<double> > > &y,
            vector<vector<double> > &y_within_mean) {
  const size_t n_chain = y.size();
  const size_t n_samp = y[0].size();
  const size_t n_var = y[0][0].size();

  // #chains x #variables
  y_within_mean = vector<vector<double> >(n_chain, vector<double>(n_var, 0.0));
  for (size_t i = 0; i < n_chain; ++i) {
    for (size_t j = 0; j < n_var; ++j) {
      for (size_t k = 0; k < n_samp; ++k) {
        y_within_mean[i][j] += y[i][k][j];
      }
      y_within_mean[i][j] /= n_samp;
    }
  }
}


void
overall_mean(const vector<vector<vector<double> > > &y,
             vector<double> &y_overall_mean) {

  const size_t n_chain = y.size();
  const size_t n_samp = y[0].size();
  const size_t n_var = y[0][0].size();
  y_overall_mean = vector<double>(n_var, 0.0); // #variables

  for (size_t i = 0; i < n_var; ++i) {
    for (size_t j = 0; j < n_chain; ++j) {
      for (size_t k = 0; k < n_samp; ++k) {
        y_overall_mean[i] += y[j][k][i];
      }
    }
    y_overall_mean[i] /= n_samp*n_chain;
  }
}


void
within_variance(const vector<vector<vector<double> > > &y,
                const vector<vector<double> > &y_within_mean,
                vector<double> &y_within_variance) {
  const size_t n_chain = y.size();
  const size_t n_samp = y[0].size();
  const size_t n_var = y[0][0].size();

  y_within_variance = vector<double>(n_var);

  for (size_t i = 0; i < n_var; ++i) {
    for (size_t j = 0; j < n_chain; ++j) {
      for (size_t k = 0; k < n_samp; ++k) {
        y_within_variance[i] += std::pow(y[j][k][i] - y_within_mean[j][i], 2);
      }
    }
    y_within_variance[i] /= n_chain*(n_samp - 1);
  }
}


void
btwn_variance(const vector<vector<double> > &y_within_mean,
              const vector<double> &y_overall_mean,
              const size_t n_samp,
              vector<double> &y_btwn_variance) {
  const size_t n_chain = y_within_mean.size();
  const size_t n_var = y_within_mean[0].size();
  y_btwn_variance = vector<double>(n_var);
  for (size_t i = 0; i < n_var; ++i) {
    for (size_t j = 0; j < n_chain; ++j) {
      y_btwn_variance[i] += std::pow(y_within_mean[j][i] - y_overall_mean[i], 2);
    }
    y_btwn_variance[i] = y_btwn_variance[i]*n_samp/(n_chain - 1);
  }
}


void
EPSR(vector<vector<mcmc_stat> > &mcmcstats,
     vector<double> &epsr) {
  epsr.clear();

  // cerr << "[in EPSR]";

  vector<vector<vector<double> > > y; // #chains x #samples x #variables
  for (size_t i = 0; i < mcmcstats.size(); ++i) { // i-th chain
    vector<vector<double> > chain;
    for (size_t j = 0; j < mcmcstats[0].size(); ++j) { // j-th sample
      vector<double> stat;
      flatten_mcmcstat(mcmcstats[i][j], stat);
      chain.push_back(stat);
    }
    y.push_back(chain);
  }

  const size_t n_samp = y[0].size();
  const size_t n_var = y[0][0].size();

  vector<vector<double> > y_within_mean; // #chain x #var
  within_mean(y, y_within_mean);

  vector<double> y_overall_mean; //#var
  overall_mean(y, y_overall_mean);

  vector<double> y_within_variance; //
  within_variance(y, y_within_mean, y_within_variance);

  vector<double> y_btwn_variance;
  btwn_variance(y_within_mean, y_overall_mean, n_samp, y_btwn_variance);

  vector<double> vhat(n_var);
  for (size_t i = 0; i < n_var; ++i) {
    vhat[i] = (y_within_variance[i]*(n_samp - 1) + y_btwn_variance[i])/n_samp;
  }

  for (size_t i = 0; i < n_var; ++i) {
    epsr.push_back(pow(vhat[i]/y_btwn_variance[i], 0.5));
  }
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////          MCMC output analysis             ///////////////////////////
///////// single chain, Batch Mean with increasing batch size //////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void
MCMC_MSE(const vector<mcmc_stat> &mcmcstats,
         const double CBM_THETA,
         const double CBM_EPS,
         double &test_val, size_t &b, size_t &a,
         bool &stop) {

  const size_t n = mcmcstats.size();
  const size_t n_sites = (mcmcstats[0].root_distr.uu +
                          mcmcstats[0].root_distr.um +
                          mcmcstats[0].root_distr.mu +
                          mcmcstats[0].root_distr.mm);

  /*determine batch number a and batch size b*/
  b = static_cast<size_t>(floor(pow(n, CBM_THETA)));
  a = static_cast<size_t>(floor(double(n)/b));

  /*whole chain average*/
  mcmc_stat ave_mcmcstats;
  vector<mcmc_stat> mcmcstats_ab(mcmcstats.end()-a*b, mcmcstats.end());
  average(mcmcstats_ab, ave_mcmcstats);

  /*batch means*/
  vector<mcmc_stat> mcmcstats_batch_means;
  for (size_t j = 0; j < a; ++j) {
    vector<mcmc_stat> batch(mcmcstats.end() - (j+1)*b, mcmcstats.end() - j*b);
    mcmc_stat batch_mean;
    average(batch, batch_mean);
    mcmcstats_batch_means.push_back(batch_mean);
  }

  /*mean squared errors*/
  mcmc_stat mse = ave_mcmcstats; // (re)initialize
  const size_t n_nodes = mcmcstats[0].start_distr.size();
  // only care about root_distr  start_distr and triad_distr
  for (size_t j = 0; j < a; ++j) {
    mse.root_distr = ((mcmcstats_batch_means[j].root_distr -
                       ave_mcmcstats.root_distr) *
                     (mcmcstats_batch_means[j].root_distr -
                      ave_mcmcstats.root_distr));
    mse.root_distr.div(double(a-1)/b);
    for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
      pair_state start_distr = ((mcmcstats_batch_means[j].start_distr[node_id] -
                                 ave_mcmcstats.start_distr[node_id]) *
                                (mcmcstats_batch_means[j].start_distr[node_id] -
                                 ave_mcmcstats.start_distr[node_id]));
      start_distr.div(double(a-1)/b);
      mse.start_distr[node_id] += start_distr;
      triple_state triad_distr = ((mcmcstats_batch_means[j].triad_distr[node_id] -
                                   ave_mcmcstats.triad_distr[node_id]) *
                                  (mcmcstats_batch_means[j].triad_distr[node_id] -
                                   ave_mcmcstats.triad_distr[node_id]));
      triad_distr.div(double(a-1)/b);
      mse.triad_distr[node_id] += triad_distr;
    }
  }

  boost::math::students_t dist(a - 1);
  double T = boost::math::quantile(complement(dist, 0.05/2));

  double cbm_max_mse = 0;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      cbm_max_mse = max(cbm_max_mse, mse.root_distr(i, j));
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        cbm_max_mse = max(cbm_max_mse, mse.root_distr(i, j));
        cbm_max_mse = max(cbm_max_mse, mse.start_distr[node_id](i, j));
        for (size_t k = 0; k < 2; ++k) {
          cbm_max_mse = max(cbm_max_mse, mse.triad_distr[node_id](i, j, k));
        }
      }
    }
  }
  test_val = pow(cbm_max_mse/(a*b), 0.5)*T;
  stop = (test_val/n_sites < CBM_EPS);
}


bool
CBM_convergence(const bool VERBOSE,
                const vector<mcmc_stat> &mcmcstats) {
  double test_val = 0;
  size_t batch_size = 0, batch_number = 0;
  bool converged;
  MCMC_MSE(mcmcstats, CBM_THETA, CBM_EPS,
           test_val, batch_size, batch_number, converged);
  if (VERBOSE)
    cerr << "\tCBM test val=" <<  test_val << "; batch_num = " << batch_number
         << "; batch_size = " << batch_size << ";";
  return converged;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////          MCMC output analysis             ///////////////////////////
///////// single chain, KL divergence                         //////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void
kl_divergence(const vector<triple_state> &P, const vector<triple_state> &Q,
              vector<double> &kld) {
  assert(P.size() == Q.size());
  kld.resize(P.size(), 0.0); // clearing not needed; values not used here
  for (size_t i = 0; i < P.size(); ++i) {
    vector<double> p, q;
    P[i].flatten(p);
    Q[i].flatten(q);
    transform(p.begin(), p.end(), p.begin(),
              bind(plus<double>(), _1, PROBABILITY_GUARD));
    transform(q.begin(), q.end(), q.begin(),
              bind(plus<double>(), _1, PROBABILITY_GUARD));
    kld[i] = kl_divergence(p, q);
  }
}

// compute KL divergence
void
kl_divergence(const mcmc_stat &P, const mcmc_stat &Q, vector<double> &kld) {
  kl_divergence(P.triad_distr, Q.triad_distr, kld);
}


// check chain convergence by KL measurement
bool
KL_convergence(const bool VERBOSE, const vector<mcmc_stat> &mcmcstats,
               const size_t keepsize, const double tol) {
  vector<double> divergence;
  vector<mcmc_stat> win1(mcmcstats.end() - 2*keepsize, mcmcstats.end() - keepsize);
  vector<mcmc_stat> win2(mcmcstats.end() - keepsize, mcmcstats.end());

  // sum statistics from two neighboring windows
  mcmc_stat sumwin1, sumwin2;
  sum(win1, sumwin1);
  sum(win2, sumwin2);

  // check how far the proportions have moved;
  kl_divergence(sumwin1, sumwin2, divergence);

  // determine convergence based on how far proportions have moved
  const double kl = *max_element(divergence.begin(), divergence.end());
 if (VERBOSE)
        cerr << "KL=" << kl << ";" ;
  return  kl < tol;
}
