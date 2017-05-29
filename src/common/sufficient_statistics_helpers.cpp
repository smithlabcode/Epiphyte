/*    Copyright (C) 2016 University of Southern California and
 *                       Andrew D. Smith and Jenny Qu
 *
 *    Authors: Jenny Qu and Andrew Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "sufficient_statistics_helpers.hpp"
#include "epiphy_utils.hpp"
#include "param_set.hpp"

#include <cassert>
#include <vector>
#include <string>

std::ostream &
operator<<(std::ostream &out, const pair_state &ps) {
  return out << ps.tostring();
}

std::ostream &
operator<<(std::ostream &out, const triple_state &ts) {
  return out << ts.tostring();
}


// compute the derivative of the transition matrix, passing values
// back through parameters
// For individual branch: T=1-exp(-branch)
static void
combine_horiz_and_vert_deriv(const param_set &ps, const size_t &node_id,
                             const pair_state &P,
                             triple_state &GP, // dim: [prev x parent x current]
                             triple_state &GP_drate, triple_state &GP_dg0,
                             triple_state &GP_dg1, triple_state &GP_dT) {
  // deriv for: (rate0, g0, g1, T)
  pair_state G(ps.g0, 1.0 - ps.g0, 1.0 - ps.g1, ps.g1);
  pair_state Q(-ps.rate0, ps.rate0, 1.0 - ps.rate0, ps.rate0 - 1.0);

  pair_state GP_denom;
  for (size_t prev = 0; prev < 2; ++prev)
    for (size_t anc = 0; anc < 2; ++anc)
      GP_denom(prev, anc) = G(prev, 0)*P(anc, 0) + G(prev, 1)*P(anc, 1);

  for (size_t prev = 0; prev < 2; ++prev)
    for (size_t anc = 0; anc < 2; ++anc)
      for (size_t cur = 0; cur < 2; ++cur)
        GP(prev, anc, cur) = G(prev, cur)*P(anc, cur)/GP_denom(prev, anc);

  pair_state dPdrate0(-ps.T[node_id], ps.T[node_id], -ps.T[node_id], ps.T[node_id]);
  pair_state dPdT(-ps.rate0, ps.rate0, 1.0 - ps.rate0, ps.rate0 - 1.0);
  pair_state dGdg0(1.0, -1.0, 0.0, 0.0);
  pair_state dGdg1(0.0, 0.0, -1.0, 1.0);

  for (size_t prev = 0; prev < 2; ++prev)
    for (size_t anc = 0; anc < 2; ++anc) {
      const double denom = GP_denom(prev, anc);
      for (size_t cur = 0; cur < 2; ++cur) {
        GP_drate(prev, anc, cur) =
          G(prev, cur)/denom*(dPdrate0(anc, cur) - P(anc, cur)/denom*
                              (G(prev, 0)*dPdrate0(anc, 0) +
                               G(prev, 1)*dPdrate0(anc, 1)));
        GP_dT(prev, anc, cur) =
          G(prev, cur)/denom*(dPdT(anc, cur) - P(anc, cur)/denom*
                              (G(prev, 0)*dPdT(anc, 0) +
                               G(prev, 1)*dPdT(anc, 1)));
        GP_dg0(prev, anc, cur) =
          (prev == 1) ? 0.0: (dGdg0(prev, cur) - G(prev, cur)/denom*
                              (dGdg0(prev, 0)*P(anc, 0) +
                               dGdg0(prev, 1)*P(anc, 1)))*P(anc, cur)/denom;
        GP_dg1(prev, anc, cur) =
          (prev == 0) ? 0.0: (dGdg1(prev, cur) - G(prev, cur)/denom*
                              (dGdg1(prev, 0)*P(anc, 0) +
                               dGdg1(prev, 1)*P(anc, 1)))*P(anc, cur)/denom;

      }
    }
}


// P(T)[0,0]: hypo -> hypo prob
// P(T)[1,1]: hyper -> hyper prob
static void
make_vertical_matrix(const double T, const double rate0, pair_state &P) {
  assert(rate0 > 0.0 && rate0 < 1.0 && T < 1.0 && T > 0.0);
  P = pair_state(1.0 - rate0*T,   rate0*T,
                 (1.0 - rate0)*T, 1.0 - (1.0 - rate0)*T);
}


void
get_transition_matrices_deriv(const param_set &ps,
                              vector<pair_state> &P,
                              vector<triple_state> &GP,
                              vector<triple_state> &GP_drate,
                              vector<triple_state> &GP_dg0,
                              vector<triple_state> &GP_dg1,
                              vector<triple_state> &GP_dT) {

  const size_t n_nodes = ps.T.size();
  P = vector<pair_state>(n_nodes);

  GP = vector<triple_state>(n_nodes); // dim: n_nodes x [prev x par x cur]
  GP_drate = vector<triple_state>(n_nodes);
  GP_dg0 = vector<triple_state>(n_nodes);
  GP_dg1 = vector<triple_state>(n_nodes);
  GP_dT = vector<triple_state>(n_nodes);

  P[0] = pair_state(1.0, 0.0, 0.0, 1.0);

  for (size_t i = 1; i < n_nodes; ++i) {
    make_vertical_matrix(ps.T[i], ps.rate0, P[i]);
    combine_horiz_and_vert_deriv(ps, i, P[i], GP[i], GP_drate[i],
                                 GP_dg0[i], GP_dg1[i], GP_dT[i]);
  }
}


static void
combine_horiz_and_vert(const pair_state &G, const pair_state &P,
                       triple_state &GP) {
  // normalization denominators
  pair_state GP_denom;
  for (size_t prev = 0; prev < 2; ++prev)
    for (size_t par = 0; par < 2; ++par)
      GP_denom(prev, par) = G(prev, 0)*P(par, 0) + G(prev, 1)*P(par, 1);

  // collecting the combined probabilities
  for (size_t prev = 0; prev < 2; ++prev)
    for (size_t par = 0; par < 2; ++par)
      for (size_t cur = 0; cur < 2; ++cur)
        GP(prev, par, cur) = G(prev, cur)*P(par, cur)/GP_denom(prev, par);
}


void
get_transition_matrices(const param_set &ps,
                        vector<pair_state> &P, vector<triple_state> &GP) {
  assert(ps.is_valid());

  const size_t n_nodes = ps.T.size();
  P = vector<pair_state>(n_nodes);
  GP = vector<triple_state>(n_nodes);
  pair_state G(ps.g0, 1.0 - ps.g0, 1.0 - ps.g1, ps.g1);

  P[0] = pair_state(1.0, 0.0, 0.0, 1.0);

  for (size_t i = 1; i < n_nodes; ++i) {
    make_vertical_matrix(ps.T[i], ps.rate0, P[i]);
    combine_horiz_and_vert(G, P[i], GP[i]);
  }
}
