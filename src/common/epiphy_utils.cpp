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

#include "epiphy_utils.hpp"
#include "MethpipeSite.hpp"
#include "PhyloTreePreorder.hpp"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <limits>

using std::vector;
using std::string;
using std::istream_iterator;
using std::pair;
using std::make_pair;
using std::numeric_limits;

static bool
parse_line(const string &line,
           vector<MSite> &sites, vector<vector<double> > &states) {

  std::istringstream iss(line);

  sites.push_back(MSite());
  if (!(iss >> sites.back().chrom >> sites.back().pos))
    return false;

  states.push_back(vector<double>(istream_iterator<double>(iss),
                                  istream_iterator<double>()));

  return true;
}

void
read_meth_table(const string &table_file,
                vector<MSite> &sites,
                vector<string> &species_names,
                vector<vector<double> > &states) {

  std::ifstream table_in(table_file.c_str());
  if (!table_in)
    throw std::runtime_error("bad table file: " + table_file);

  string line;
  getline(table_in, line);
  std::istringstream iss(line);
  copy(istream_iterator<string>(iss), istream_iterator<string>(),
       std::back_inserter(species_names));

  while (getline(table_in, line))
    if (line[0] != '#')
      if (!parse_line(line, sites, states) ||
          states.back().size() != species_names.size())
        throw std::runtime_error("bad table file line: " + line);
}


void
mark_useable_sites(const vector<size_t> &subtree_sizes,
                   const vector<MSite> &sites,
                   const vector<vector<double> > &tree_probs,
                   const size_t desert_size, const size_t node_id,
                   vector<vector<bool> > &marks) {

  if (is_leaf(subtree_sizes[node_id])) {
    size_t prev_obs = numeric_limits<size_t>::max();

    for (size_t i = 0; i < tree_probs.size();) {
      // find next observed site
      while(i < tree_probs.size() &&
            missing_meth_value(tree_probs[i][node_id])) ++i;

      if (i < tree_probs.size()) {
        marks[i][node_id] = true;

        // fill in interval if not end of a desert
        if (i < tree_probs.size() && prev_obs < i &&
            distance(sites[prev_obs], sites[i]) <= desert_size)
          for (size_t j = prev_obs + 1; j < i; ++j)
            marks[j][node_id] = true;

        prev_obs = i;
        ++i;
      }
    }
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      mark_useable_sites(subtree_sizes, sites, tree_probs,
                         desert_size, node_id + count, marks);
      // non-desert sites in an internal node
      // is the intersection of non-deserts of all its children
      if (count == 1) {
        for (size_t i = 0; i < tree_probs.size(); ++i)
          marks[i][node_id] = marks[i][node_id + count];
      } else {
        for (size_t i = 0; i < tree_probs.size(); ++i)
          marks[i][node_id] = marks[i][node_id] && marks[i][node_id + count];
      }
      count += subtree_sizes[node_id + count];
    }
  }
}

void
mark_useable_sites(const vector<size_t> subtree_sizes,
                   const vector<MSite> &sites,
                   const vector<vector<double> > &tree_probs,
                   const size_t desert_size,
           vector<vector<bool> > &marks) {
  marks = vector<vector<bool> >(tree_probs.size(),
                                vector<bool>(subtree_sizes.size(), false));
  mark_useable_sites(subtree_sizes, sites, tree_probs, desert_size, 0, marks);
}


void
marks_to_blocks(const vector<size_t> subtree_sizes,
                const vector<MSite> &sites,
                const size_t desert_size,
                const vector<vector<bool> > &marks,
                vector<vector<pair<size_t, size_t> > > &blocks) {
  blocks = vector<vector<pair<size_t, size_t> > >(subtree_sizes.size());
  for (size_t node_id = 0; node_id < subtree_sizes.size(); ++node_id) {
    for (size_t i = 0; i < marks.size(); ++i) {
      if ((i == 0 && marks[i][node_id]) ||
          (i > 0 && !marks[i - 1][node_id] && marks[i][node_id]) ||
          (i > 0 && marks[i - 1][node_id] && marks[i][node_id] &&
           distance(sites[i - 1], sites[i]) > desert_size)) {
        blocks[node_id].push_back(std::make_pair(i, i));
      } else if (i > 0 && marks[i - 1][node_id] && marks[i][node_id]) {
        blocks[node_id].back().second = i;
      }
    }
  }
}
