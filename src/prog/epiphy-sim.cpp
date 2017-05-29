/*****************************************************************************
 *  epiphy-sim: a program to simulate methylation states for species in
 *  a phylogenetic tree according to inheritance and auto-correlation
 *  of methylation states, and WGBS data sets for extant species.
 *
 *  Copyright (C) 2016 University of Southern California and
 *                     Jenny Qu and Andrew D Smith
 *
 *  Authors: Jenny Qu and Andrew D Smith
 *
 *  This program is free software: you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see
 *  <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <sstream>
#include <random>
#include <iomanip>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "MethpipeSite.hpp"

#include "PhyloTreePreorder.hpp"
#include "param_set.hpp"
#include "sufficient_statistics_helpers.hpp"
#include "optimize_params.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::ostream_iterator;

using std::max;
using std::min;
using std::pair;

/*
  f0: u->u transition probability in G
  f1: m->m transition probability in G
  g0: u->u transition probability in GP
  g1: m->m transition probability in GP
  pi0: initial probability of u
  states: u=0 m=1
 */

static void
separate_regions(const size_t desert_size, vector<MSite> &sites,
                 vector<pair<size_t, size_t> > &blocks) {
  for (size_t i = 0; i < sites.size(); ++i)
    if (i == 0 || distance(sites[i - 1], sites[i]) > desert_size)
      blocks.push_back(std::make_pair(i, i));
    else blocks.back().second = i;
}


class single_edge_sampler {
public:
  single_edge_sampler(const pair_state &P) { // need to allow seed
    dis = std::uniform_real_distribution<>(0, 1);
    std::random_device rd;
    gen = std::mt19937_64(rd());

    // JQU: why? Can P(0, 0) + P(0, 1) not sum to 1?
    m_gvn_parent_u = P(0, 1)/max(1.0, P(0, 0) + P(0, 1));
    m_gvn_parent_m = P(1, 1)/max(1.0, P(1, 0) + P(1, 1));
  }
  bool operator()(const bool parent_is_m) const {
    return dis(gen) < (parent_is_m ? m_gvn_parent_m : m_gvn_parent_u);
  }
  string tostring() const {
    std::ostringstream oss;
    oss << std::setprecision(3);
    // organized set to show the "G=(g0, g1)" order
    oss << '['
        << 1.0 - m_gvn_parent_u << ',' << m_gvn_parent_u << "]["
        << 1.0 - m_gvn_parent_m << ',' << m_gvn_parent_m << ']';
    return oss.str();
  }
private:
  mutable std::uniform_real_distribution<> dis;
  mutable std::mt19937_64 gen;
  double m_gvn_parent_u;
  double m_gvn_parent_m;
};

std::ostream&
operator<<(std::ostream &out, const single_edge_sampler &ses) {
  return out << ses.tostring();
}


class two_edge_sampler {
public:
  two_edge_sampler(const triple_state &GP) { // need to allow seed
    dis = std::uniform_real_distribution<>(0, 1);
    std::random_device rd;
    gen = std::mt19937_64(rd());

    m_gvn_left_u_par_u = GP(0, 0, 1)/max(1.0, GP(0, 0, 0) + GP(0, 0, 1));
    m_gvn_left_u_par_m = GP(0, 1, 1)/max(1.0, GP(0, 1, 0) + GP(0, 1, 1));
    m_gvn_left_m_par_u = GP(1, 0, 1)/max(1.0, GP(1, 0, 0) + GP(1, 0, 1));
    m_gvn_left_m_par_m = GP(1, 1, 1)/max(1.0, GP(1, 1, 0) + GP(1, 1, 1));
  }
  bool operator()(const bool left_is_m, const bool par_is_m) const {
    // the ">" below is because the "true" value is for methylated
    return dis(gen) < (left_is_m ?
                       (par_is_m ? m_gvn_left_m_par_m : m_gvn_left_m_par_u) :
                       (par_is_m ? m_gvn_left_u_par_m : m_gvn_left_u_par_u));
  }
  string tostring() const {
    std::ostringstream oss;
    oss << std::setprecision(3);
    oss << '['
        << 1.0 - m_gvn_left_u_par_u << ',' << m_gvn_left_u_par_u << "]["
        << 1.0 - m_gvn_left_u_par_m << ',' << m_gvn_left_u_par_m << "]\n["
        << 1.0 - m_gvn_left_m_par_u << ',' << m_gvn_left_m_par_u << "]["
        << 1.0 - m_gvn_left_m_par_m << ',' << m_gvn_left_m_par_m << ']';
    return oss.str();
  }
private:
  mutable std::uniform_real_distribution<> dis;
  mutable std::mt19937_64 gen;

  double m_gvn_left_u_par_u;
  double m_gvn_left_u_par_m;
  double m_gvn_left_m_par_u;
  double m_gvn_left_m_par_m;
};

std::ostream&
operator<<(std::ostream &out, const two_edge_sampler &tes) {
  return out << tes.tostring();
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* given state at node[node_id], simulate nodes in the subtree */
template <class T>
static void
simulate_site_start(const vector<size_t> &subtree_sizes,
                    const vector<single_edge_sampler> &Psamp,
                    const size_t node_id,
                    vector<pair_state> &start_counts, vector<T> &states) {
  if (!is_leaf(subtree_sizes[node_id])) {
    const T current_state = states[node_id];
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      states[child_id] = Psamp[child_id](current_state);
      start_counts[child_id](current_state, states[child_id])++;
      simulate_site_start(subtree_sizes, Psamp, child_id, start_counts, states);
      count += subtree_sizes[child_id];
    }
  }
}

/* simulate root node, then recursively simulate subtrees */
template <class T>
static void
simulate_site_start(const vector<size_t> &subtree_sizes, const param_set &ps,
                    const vector<single_edge_sampler> &Psamp,
                    pair<double, double> &root_start_counts,
                    vector<pair_state> &start_counts,
                    vector<T> &states) {

  states.resize(subtree_sizes.size());

  // sample root at start position
  std::random_device rd;
  std::mt19937_64 gen(rd());
  states[0] = std::uniform_real_distribution<>(0, 1)(gen) > ps.pi0; // == [< pi1]

  root_start_counts.first += !states[0]; // unmeth state count
  root_start_counts.second += states[0];

  simulate_site_start(subtree_sizes, Psamp, 0, start_counts, states);
}

template <class T>
static void
simulate_site(const vector<size_t> &subtree_sizes,
              const vector<two_edge_sampler> &GPsamp,
              const vector<bool> &prev_states, const size_t node_id,
              vector<triple_state> &triad_counts,
              vector<T> &states) {

  if (!is_leaf(subtree_sizes[node_id])) {
    const T current_state = states[node_id];
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      states[child_id] = GPsamp[child_id](prev_states[child_id], current_state);
      triad_counts[child_id](prev_states[child_id],
                             current_state, states[child_id])++;
      simulate_site(subtree_sizes, GPsamp, prev_states, child_id,
                    triad_counts, states);
      count += subtree_sizes[child_id];
    }
  }
}

template <class T>
static void
simulate_site(const vector<size_t> &subtree_sizes,
              const single_edge_sampler &Gsamp,
              const vector<two_edge_sampler> &GPsamp,
              const vector<bool> &prev_states,
              pair_state &root_counts,
              vector<triple_state> &triad_counts,
              vector<T> &states) {

  states.resize(subtree_sizes.size());
  states[0] = Gsamp(prev_states[0]); // sample root using state to left

  root_counts(prev_states[0], states[0])++;

  simulate_site(subtree_sizes, GPsamp, prev_states, 0, triad_counts, states);
}


int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;
    string leaf_only_file;
    size_t desert_size = 1000;
    bool independent_sites = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "simulate methylomes "
                           "related by a phylogenetic tree",
                           "<sites-file> <params-file>");
    opt_parse.add_opt("indep", 'I', "sites are independent",
                      false, independent_sites);
    opt_parse.add_opt("desert", 'd',
                      "desert size (default " + std::to_string(1000) + ")",
                      false, desert_size);
    opt_parse.add_opt("leaf-file", 'L',
                      "write leaf-only data to this file",
                      false, leaf_only_file);
    opt_parse.add_opt("output", 'o', "name of output file "
                      "(default: stdout)", true, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file(leftover_args.front());
    const string param_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[reading parameters]" << endl;
    PhyloTreePreorder t;
    param_set ps;
    ps.read(param_file, t);
    if (VERBOSE)
      cerr << "==> parameters={" << ps << "}" << endl;

    vector<size_t> subtree_sizes;
    t.get_subtree_sizes(subtree_sizes);
    vector<string> leaf_names;
    t.get_leaf_names(leaf_names);
    t.assign_missing_node_names();

    const size_t n_nodes = t.get_size();

    if (VERBOSE)
      cerr << "[computing transition matrices]" << endl;

    pair_state F(ps.f0, 1.0 - ps.f0, 1.0 - ps.f1, ps.f1);
    vector<pair_state> P;
    vector<triple_state> GP;
    get_transition_matrices(ps, P, GP);

    single_edge_sampler Gsamp(F);
    vector<single_edge_sampler> Psamp;
    vector<two_edge_sampler> GPsamp;
    for (size_t i = 0; i < n_nodes; ++i) {
      Psamp.push_back(single_edge_sampler(P[i]));
      GPsamp.push_back(two_edge_sampler(GP[i]));
    }

    if (VERBOSE) {
      cerr << "==> horizontal transitions:" << endl
           << Gsamp << endl;
      cerr << "==> vertical transitions:" << endl;
      for (size_t i = 0; i < n_nodes; ++i)
        cerr << Psamp[i] << endl;
      cerr << "==> combined transitions:" << endl;
      for (size_t i = 0; i < n_nodes; ++i)
        cerr << GPsamp[i] << endl;
    }

    if (VERBOSE)
      cerr << "[loading sites]" << endl;
    std::ifstream sites_in(cpgs_file.c_str());
    if (!sites_in)
      throw std::runtime_error("bad input file: " + cpgs_file);
    vector<MSite> sites;
    MSite site;
    while (sites_in >> site) sites.push_back(site);
    const size_t n_sites = sites.size();
    if (VERBOSE)
      cerr << "==> n_sites=" << n_sites << endl;

    pair<double, double> root_start_counts;
    pair_state root_counts;
    vector<pair_state> start_counts(n_nodes);
    vector<triple_state> triad_counts(n_nodes);

    vector<vector<bool> > states(n_sites);

    if (independent_sites) {
      if (VERBOSE)
        cerr << "[simulating]" << endl;
      for (size_t i = 0; i < states.size(); ++i)
        simulate_site_start(subtree_sizes, ps, Psamp,
                            root_start_counts, start_counts, states[i]);
    }
    else {
      if (VERBOSE)
        cerr << "[separating by deserts]" << endl;
      vector<pair<size_t, size_t> > blocks;
      separate_regions(desert_size, sites, blocks);
      if (VERBOSE)
        cerr << "==> n_resets=" << blocks.size() << endl;

      if (VERBOSE)
        cerr << "[simulating]" << endl;

      for (size_t i = 0; i < blocks.size(); ++i) {

        const size_t start = blocks[i].first;
        const size_t end = blocks[i].second;

        simulate_site_start(subtree_sizes, ps, Psamp,
                            root_start_counts, start_counts, states[start]);
        for (size_t pos = start + 1; pos <= end; ++pos)
          simulate_site(subtree_sizes, Gsamp, GPsamp, states[pos - 1],
                        root_counts, triad_counts, states[pos]);
      }
    }

    vector<size_t> leaves_preorder;
    subtree_sizes_to_leaves_preorder(subtree_sizes, leaves_preorder);

    std::ofstream out(outfile.c_str());
    vector<string> node_names;
    t.get_node_names(node_names);
    copy(node_names.begin(), node_names.end(),
         ostream_iterator<string>(out, "\t"));
    out << endl;
    for (size_t i = 0; i < states.size(); ++i) {
      out << sites[i].chrom << '\t' << sites[i].pos;
      for (size_t j = 0; j < states[i].size(); ++j)
        out << '\t' << states[i][j];
      out << endl;
    }

    if (VERBOSE) {
      cerr << "root_start_counts:\n"
           << root_start_counts.first/(root_start_counts.first +
                                       root_start_counts.second) << endl
           << "root_counts:\n" << root_counts << endl
           << "start_counts:\n";
      copy(start_counts.begin(), start_counts.end(),
           ostream_iterator<pair_state>(cerr, "\n"));
      if (!independent_sites) {
        cerr << "triad_counts:\n";
        copy(triad_counts.begin(), triad_counts.end(),
             ostream_iterator<triple_state>(cerr, "\n"));
      }
      const double llk = independent_sites ?
        log_likelihood(subtree_sizes, ps, root_start_counts,
                       root_counts, start_counts) :
        log_likelihood(subtree_sizes, ps, root_start_counts,
                       root_counts, start_counts, triad_counts);
      cerr << "log_likelihood=" << llk << endl;
    }

    if (!leaf_only_file.empty()) {
      std::ofstream leaf_out(leaf_only_file.c_str());

      vector<string> leaf_names;
      t.get_leaf_names(leaf_names);
      copy(leaf_names.begin(), leaf_names.end(),
           ostream_iterator<string>(leaf_out, "\t"));
      leaf_out << endl;

      for (size_t i = 0; i < states.size(); ++i) {
        leaf_out << sites[i].chrom << '\t' << sites[i].pos;
        for (size_t j = 0; j < states[i].size(); ++j)
          if (is_leaf(subtree_sizes[j]))
            leaf_out << '\t' << states[i][j];
        leaf_out << endl;
      }
    }

  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
