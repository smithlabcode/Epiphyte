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

#ifndef EPIPHY_UTILS_HPP
#define EPIPHY_UTILS_HPP

#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <numeric>
#include <iterator>

template <class T> double
log_sum_log(const T p, const T q) {

  if (p == 0) return q;
  if (q == 0) return p;

  const T larger = std::max(p, q);
  return larger + std::log(1.0 + exp(std::min(p, q) - larger));
}

template <class T> double
kronecker_delta(const T &a, const T &b) {
  return (a == b) ? 1.0 : 0.0;
}

inline double
sign(const double val) {
  return (val == 0.0) ? 0.0 : (val > 0.0 ? 1.0 : -1.0);
}

inline double
kl_divergence(const std::vector<double> &P, const std::vector<double> &Q) {
  assert(P.size()==Q.size());
  const double psum = std::accumulate(P.begin(), P.end(), 0.0);
  const double qsum = std::accumulate(Q.begin(), Q.end(), 0.0);
  const double norm = std::log(psum) - std::log(qsum);

  double d = 0.0;
  for (size_t i = 0; i < P.size(); ++i)
    d += P[i]/psum*(std::log(P[i]) - std::log(Q[i]) - norm);
  return d;
}

inline bool
missing_meth_value(const double x) {return x == -1.0;}

inline bool
valid_probability(const double x) {return x >= 0.0 && x <= 1.0;}

inline bool
valid_meth_value(const double x) {
  return valid_probability(x) || missing_meth_value(x);
}

inline double
away_from_extremes(const double x, const double epsilon) {
  return std::min(std::max(epsilon, x), 1.0 - epsilon);
}

class MSite;

void
read_meth_table(const std::string &table_file,
                std::vector<MSite> &sites,
                std::vector<std::string> &species_names,
                std::vector<std::vector<double> > &states);

void
mark_useable_sites(const std::vector<size_t> subtree_sizes,
                   const std::vector<MSite> &sites,
                   const std::vector<std::vector<double> > &tree_probs,
                   const size_t desert_size,
                   std::vector<std::vector<bool> > &marks);


#endif
