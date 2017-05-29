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

#ifndef SUFFICIENT_STATISTICS_HELPERS_HPP
#define SUFFICIENT_STATISTICS_HELPERS_HPP

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////  STRUCTS FOR HOLDING PAIRS AND TRIPLES OF STATES  ////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include <MethpipeSite.hpp>

#include <vector>
#include <string>
#include <cmath>
#include <sstream>

struct param_set;

struct pair_state {
  double uu, um; // think 2x2 matrix
  double mu, mm;

  pair_state(double _uu, double _um, double _mu, double _mm) :
    uu(_uu), um(_um), mu(_mu), mm(_mm) {}
  pair_state() : uu(0.0), um(0.0), mu(0.0), mm(0.0) {}

  // WARNING: no range checking (for i,j > 1)
  double operator()(int i, int j) const { // version for l-values
    return (i == 0) ? (j == 0 ? uu : um) : (j == 0 ? mu : mm);
  }
  double & operator()(int i, int j) { // version for r-values
    return (i == 0) ? (j == 0 ? uu : um) : (j == 0 ? mu : mm);
  }
  pair_state operator+(const pair_state &other) const {
    return pair_state(uu + other.uu, um + other.um,
                      mu + other.mu, mm + other.mm);
  }
  pair_state operator-(const pair_state &other) const {
    return pair_state(uu - other.uu, um - other.um,
                      mu - other.mu, mm - other.mm);
  }
  pair_state operator*(const pair_state &other) const {
    return pair_state(uu*other.uu, um*other.um,
                      mu*other.mu, mm*other.mm);
  }
  void operator+=(const pair_state &other) {
    uu += other.uu; um += other.um;
    mu += other.mu; mm += other.mm;
  }
  void operator/=(const pair_state &other) {
    uu /= other.uu; um /= other.um;
    mu /= other.mu; mm /= other.mm;
  }
  void div(const double x) {
    uu /= x; um /= x;
    mu /= x; mm /= x;
  }
  void to_probabilities() {
    const double u_denom = uu + um;
    uu /= u_denom;
    um /= u_denom;
    const double m_denom = mu + mm;
    mu /= m_denom;
    mm /= m_denom;
  }
  void make_logs() {
    uu = std::log(uu); um = std::log(um);
    mu = std::log(mu); mm = std::log(mm);
  }
  void flatten(std::vector<double> &p) const {
    p.clear();
    p.push_back(uu); p.push_back(um);
    p.push_back(mu); p.push_back(mm);
  }
  std::string tostring() const {
    std::ostringstream oss;
    oss << "[" << uu << ", " << um << "]\n"
        << "[" << mu << ", " << mm << "]";
    return oss.str();
  }
};

std::ostream &
operator<<(std::ostream &out, const pair_state &ps);

struct triple_state {
  double uuu, uum; // think: a 2x2 matrix
  double umu, umm;

  double muu, mum; // think: another 2x2 matrix
  double mmu, mmm;

  triple_state() : uuu(0.0), uum(0.0), umu(0.0), umm(0.0),
                   muu(0.0), mum(0.0), mmu(0.0), mmm(0.0) {};
  triple_state(double _uuu, double _uum, double _umu, double _umm,
               double _muu, double _mum, double _mmu, double _mmm) :
    uuu(_uuu), uum(_uum), umu(_umu), umm(_umm),
    muu(_muu), mum(_mum), mmu(_mmu), mmm(_mmm) {};

  // WARNING: no range checking
  double operator()(int i, int j, int k) const { // version for r-value (rhs)
    return  i == 0 ?
      (j == 0 ? (k == 0 ? uuu : uum) : (k == 0 ? umu : umm)) :
      (j == 0 ? (k == 0 ? muu : mum) : (k == 0 ? mmu : mmm));
  }
  double & operator()(int i, int j, int k) { // version for l-value (lhs)
    return  i == 0 ?
      (j == 0 ? (k == 0 ? uuu : uum) : (k == 0 ? umu : umm)) :
      (j == 0 ? (k == 0 ? muu : mum) : (k == 0 ? mmu : mmm));
  }
  triple_state operator+(const triple_state &other) const {
    return triple_state(uuu + other.uuu, uum + other.uum,
                        umu + other.umu, umm + other.umm,
                        muu + other.muu, mum + other.mum,
                        mmu + other.mmu, mmm + other.mmm);
  }
  triple_state operator-(const triple_state &other) const {
    return triple_state(uuu - other.uuu, uum - other.uum,
                        umu - other.umu, umm - other.umm,
                        muu - other.muu, mum - other.mum,
                        mmu - other.mmu, mmm - other.mmm);
  }
  triple_state operator*(const triple_state &other) const {
    return triple_state(uuu*other.uuu, uum*other.uum,
                        umu*other.umu, umm*other.umm,
                        muu*other.muu, mum*other.mum,
                        mmu*other.mmu, mmm*other.mmm);
  }
  void operator+=(const triple_state &other) {
    uuu += other.uuu; uum += other.uum;
    umu += other.umu; umm += other.umm;
    muu += other.muu; mum += other.mum;
    mmu += other.mmu; mmm += other.mmm;
  }
  void operator/=(const triple_state &other) {
    uuu /= other.uuu; uum /= other.uum;
    umu /= other.umu; umm /= other.umm;
    muu /= other.muu; mum /= other.mum;
    mmu /= other.mmu; mmm /= other.mmm;
  }
  void to_probabilities() {
    const double uu_denom = uuu + uum;
    uuu /= uu_denom;
    uum /= uu_denom;
    const double um_denom = umu + umm;
    umu /= um_denom;
    umm /= um_denom;
    const double mu_denom = muu + mum;
    muu /= mu_denom;
    mum /= mu_denom;
    const double mm_denom = mmu + mmm;
    mmu /= mm_denom;
    mmm /= mm_denom;
  }
  void div(const double x) {
    uuu /= x; uum /= x;
    umu /= x; umm /= x;

    muu /= x; mum /= x;
    mmu /= x; mmm /= x;
  }
  void make_logs() {
    uuu = std::log(uuu); uum = std::log(uum);
    umu = std::log(umu); umm = std::log(umm);

    muu = std::log(muu); mum = std::log(mum);
    mmu = std::log(mmu); mmm = std::log(mmm);
  }
  void flatten(std::vector<double> &p) const {
    p.clear();
    p.push_back(uuu); p.push_back(uum);
    p.push_back(umu); p.push_back(umm);
    p.push_back(muu); p.push_back(mum);
    p.push_back(mmu); p.push_back(mmm);
  }
  std::string tostring() const {
    std::ostringstream oss;
    oss << "u [" << uuu << ", " << uum << "]\n"
        << "  [" << umu << ", " << umm << "]\n"
        << "m [" << muu << ", " << mum << "]\n"
        << "  [" << mmu << ", " << mmm << "]";
    return oss.str();
  }
};

std::ostream &
operator<<(std::ostream &out, const triple_state &ts);

void
get_transition_matrices(const param_set &ps,
                        std::vector<pair_state> &P,
                        std::vector<triple_state> &GP);

void
get_transition_matrices_deriv(const param_set &ps,
                              std::vector<pair_state> &P,
                              std::vector<triple_state> &GP,
                              std::vector<triple_state> &GP_drate,
                              std::vector<triple_state> &GP_dg0,
                              std::vector<triple_state> &GP_dg1,
                              std::vector<triple_state> &GP_dT);

/* The count_triads function resides in the header and not cpp file
   because it is templated. */
template <class T>
static void
count_triads(const std::vector<size_t> &subtree_sizes,
             const std::vector<size_t> &parent_ids,
             const std::vector<std::vector<T> > &tree_states,
             const std::vector<std::pair<size_t, size_t> > &reset_points,
             std::pair<double, double> &root_start_counts,
             pair_state &root_counts,
             std::vector<pair_state> &start_counts,
             std::vector<triple_state> &triad_counts) {

  root_start_counts = std::make_pair(0.0, 0.0);
  triad_counts = std::vector<triple_state>(subtree_sizes.size());
  start_counts = std::vector<pair_state>(subtree_sizes.size());
  root_counts = pair_state();

  for (size_t i = 0; i < reset_points.size(); ++i) {

    const size_t start = reset_points[i].first;
    const size_t end = reset_points[i].second;

    root_start_counts.first += !tree_states[start][0];
    root_start_counts.second += tree_states[start][0];

    for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
      const size_t parent = tree_states[start][parent_ids[node_id]];
      const size_t curr = tree_states[start][node_id];
      start_counts[node_id](parent, curr) += 1.0;
    }

    for (size_t pos = start + 1; pos <= end; ++pos) {

      root_counts(tree_states[pos - 1][0], tree_states[pos][0])++;

      for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
        const size_t parent = tree_states[pos][parent_ids[node_id]];
        const size_t prev = tree_states[pos - 1][node_id];
        const size_t curr = tree_states[pos][node_id];
        triad_counts[node_id](prev, parent, curr) += 1.0;
      }
    }
  }
}


// general case (marks)
template <class T>
static void
count_triads(const std::vector<size_t> &subtree_sizes,
             const std::vector<size_t> &parent_ids,
             const std::vector<std::vector<T> > &tree_states,
             const std::vector<std::vector<bool> > &marks,
             const std::vector<MSite> &sites,
             const size_t desert_size,
             std::pair<double, double> &root_start_counts,
             pair_state &root_counts,
             std::vector<pair_state> &start_counts,
             std::vector<triple_state> &triad_counts) {

  const size_t n_sites = tree_states.size();
  const size_t n_nodes = subtree_sizes.size();

  root_start_counts = std::make_pair(0.0, 0.0);
  triad_counts = std::vector<triple_state>(n_nodes);
  start_counts = std::vector<pair_state>(n_nodes);
  root_counts = pair_state();

  for (size_t i = 0; i < n_sites; ++i) {
    for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
      if (marks[i][node_id]) {
        bool has_prev = (i > 0 && marks[i-1][node_id] &&
                         distance(sites[i-1], sites[i]) <= desert_size);
        bool has_par = (node_id > 0 &&  marks[i][parent_ids[node_id]]);

        if (has_prev) { // has prev
          if (has_par) { // has parent
            const size_t parent = tree_states[i][parent_ids[node_id]];
            const size_t prev = tree_states[i - 1][node_id];
            const size_t curr = tree_states[i][node_id];
            triad_counts[node_id](prev, parent, curr) += 1.0;
          } else { //no parent
            const size_t prev = tree_states[i - 1][node_id];
            const size_t curr = tree_states[i][node_id];
            root_counts(prev, curr) += 1.0;
          }
        } else { // no prev
          if (has_par) { // has parent
            const size_t parent = tree_states[i][parent_ids[node_id]];
            const size_t curr = tree_states[i][node_id];
            start_counts[node_id](parent, curr) += 1.0;
          } else { //no parent
            root_start_counts.first += !tree_states[i][node_id];
            root_start_counts.second += tree_states[i][node_id];
          }
        }
      }
    }
  }

  // // tot number of different types of counts/sites
  // size_t root_start_tot;
  // size_t root_tot;
  // std::vector<size_t> start_tot;
  // std::vector<size_t> triad_tot;

  // root_tot = (root_counts(0, 0) + root_counts(0, 1) +
  //             root_counts(1, 0) + root_counts(1, 1));
  // root_start_tot = root_start_counts.first + root_start_counts.second;
  // start_tot = vector<size_t>(subtree_sizes.size(), 0);
  // triad_tot = vector<size_t>(subtree_sizes.size(), 0);
  // for (size_t node_id = 1; node < n_nodes; ++node_id) {
  //   tot_triad[node_id] =
  //     (triad_counts[node_id](0, 0, 0) + triad_counts[node_id](0, 0, 1) +
  //      triad_counts[node_id](0, 1, 0) + triad_counts[node_id](0, 1, 1) +
  //      triad_counts[node_id](1, 0, 0) + triad_counts[node_id](1, 0, 1) +
  //      triad_counts[node_id](1, 1, 0) + triad_counts[node_id](1, 1, 1));
  //   start_tot[node_id] =
  //     (start_counts[node_id](0, 0) + start_counts[node_id](0, 0) +
  //      start_counts[node_id](0, 1) + start_counts[node_id](1, 1));
  // }
}

#endif
