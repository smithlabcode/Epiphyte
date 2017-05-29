/*    Copyright (C) 2015-16 University of Southern California and
 *                          Andrew D. Smith and Jenny Qu
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

/******************************************************************************/

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <random>
#include <algorithm>  //std::max, min, random_shuffle
#include <cmath>      //std::abs
#include <limits>     //std::numeric_limits
#include <iterator>   //std::distance
#include <random>

/* from smithlab_cpp */
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

/* from methpipe */
#include "MethpipeSite.hpp"

/* headers for epigenomic evolution */
#include "PhyloTreePreorder.hpp"
#include "PhyloTree.hpp"
#include "param_set.hpp"
#include "epiphy_utils.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::numeric_limits;
using std::min;
using std::max;
using std::istringstream;
using std::to_string;
using std::ostream_iterator;


static void
separate_regions(const size_t desert_size, const vector<MSite> &sites,
                 vector<pair<size_t, size_t> > &blocks) {
  for (size_t i = 0; i < sites.size(); ++i)
    if (i == 0 || distance(sites[i - 1], sites[i]) > desert_size)
      blocks.push_back(std::make_pair(i, i));
    else blocks.back().second = i;
}


static void
marks_to_blocks(const size_t desert_size,
                const vector<size_t> &subtree_sizes,
                const vector<MSite> &sites,
                const vector<vector<bool> > &marks,
                vector<vector<pair<size_t, size_t> > > &blocks_by_sp) {
  const size_t n_nodes = subtree_sizes.size();
  const size_t n_sites = sites.size();
  blocks_by_sp = vector<vector<pair<size_t, size_t> > >(n_nodes);
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    for (size_t i = 0; i < n_sites; ++i) {
      if (marks[i][node_id]) {
        if (!blocks_by_sp[node_id].size()) {
          blocks_by_sp[node_id].push_back(std::make_pair(i, i));
        } else {
          const size_t last_site = blocks_by_sp[node_id].back().second;
          if (distance(sites[last_site], sites[i]) > desert_size)
            blocks_by_sp[node_id].push_back(std::make_pair(i, i));
          else
            blocks_by_sp[node_id].back().second = i;
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////   Build domains for each species    /////////////////////
////////////////////////////////////////////////////////////////////////////////
// keep foreground domain scores (sum of 1-methprob)
static void
get_domain_scores(const vector<size_t> &subtree_sizes,
                  const vector<vector<double> > &meth_probs,
                  const vector<pair<size_t, size_t> > &blocks,
                  vector<vector<double> > &domain_scores) {
  const size_t n_nodes = subtree_sizes.size();
  vector<vector<bool> > classes(n_nodes);

  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    for (size_t i = 0; i < blocks.size(); ++i) {
      const size_t start = blocks[i].first;
      const size_t end = blocks[i].second;
      for (size_t pos = start; pos <= end; ++pos) {
        double p = meth_probs[pos][node_id];
        bool state = p > 0.5;
        if (pos == start) {
          classes[node_id].push_back(state);
          if (!state)
            domain_scores[node_id].push_back(0.5 - p);
        } else {
          if (classes[node_id].back() == state) {
            if (!state)
              domain_scores[node_id].back() += 0.5 - p;
          } else {
            classes[node_id].push_back(state);
            if (!state)
              domain_scores[node_id].push_back(0.5 - p);
          }
        }
      }
    }
  }
}


// general case with marked sites
// (blocks defined separately for each node)
static void
get_domain_scores(const vector<size_t> &subtree_sizes,
                  const vector<vector<double> > &meth_probs,
                  const vector<vector<pair<size_t, size_t> > > &blocks_by_sp,
                  vector<vector<double> > &domain_scores) {

  const size_t n_nodes = subtree_sizes.size();
  vector<vector<bool> > classes(n_nodes);

  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    for (size_t i = 0; i < blocks_by_sp[node_id].size(); ++i) {
      const size_t start = blocks_by_sp[node_id][i].first;
      const size_t end = blocks_by_sp[node_id][i].second;
      for (size_t pos = start; pos <= end; ++pos) {
        double p = meth_probs[pos][node_id];
        bool state = p > 0.5;
        if (pos == start) {
          classes[node_id].push_back(state);
          if (!state)
            domain_scores[node_id].push_back(0.5 - p);
        } else {
          if (classes[node_id].back() == state) {
            if (!state)
              domain_scores[node_id].back() += 0.5 - p;
          } else {
            classes[node_id].push_back(state);
            if (!state)
              domain_scores[node_id].push_back(0.5 - p);
          }
        }
      }
    }
  }
}


static void
shuffle_cpgs(const size_t seed,
             const vector<size_t> &subtree_sizes,
             const vector<pair<size_t, size_t> > &blocks,
             vector<vector<double> > meth_probs,
             vector<double> &random_scores) {
  const size_t n_nodes = subtree_sizes.size();
  vector<vector<double> > domain_scores(n_nodes);
  srand(seed);
  random_shuffle(meth_probs.begin(), meth_probs.end());
  vector<double> scores;
  get_domain_scores(subtree_sizes, meth_probs, blocks, domain_scores);

  bool DEBUG = false;
  if (DEBUG) {
    string outfile = "shuffled_scores_tmp.txt";
    std::ofstream of;
    of.open(outfile.c_str());
    std::ostream out(of.rdbuf());
    for (size_t node_id = 0; node_id < n_nodes; ++node_id)
      for (size_t i = 0; i < domain_scores[node_id].size(); ++i)
        out << node_id << "\t" << domain_scores[node_id][i] << endl;
  }

  for (size_t node_id = 0; node_id < n_nodes; ++node_id)
    random_scores.insert(random_scores.end(), domain_scores[node_id].begin(),
                         domain_scores[node_id].end());
  sort(random_scores.begin(), random_scores.end());
}


//shuffle within each block in each species
static void
shuffle_cpgs(const size_t seed,
             const vector<size_t> &subtree_sizes,
             const vector<vector<pair<size_t, size_t> > > &blocks_by_sp,
             vector<vector<double> > meth_probs,
             vector<double> &random_scores) {
  const size_t n_nodes = subtree_sizes.size();
  const size_t n_sites = meth_probs.size();
  vector<vector<double> > domain_scores(n_nodes);
  srand(seed);

  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    vector<double> nonmissing_meth;
    for (size_t i = 0; i < n_sites; ++i) {
      if (!missing_meth_value(meth_probs[i][node_id]))
        nonmissing_meth.push_back(meth_probs[i][node_id]);
    }
    random_shuffle(nonmissing_meth.begin(), nonmissing_meth.end());
    size_t pos = 0;
    size_t count = 0;
    while (pos < n_sites) {
      if (!missing_meth_value(meth_probs[pos][node_id])) {
        meth_probs[pos][node_id] = nonmissing_meth[count];
        ++count;
      }
      ++pos;
    }
  }
  vector<double> scores;
  get_domain_scores(subtree_sizes, meth_probs, blocks_by_sp, domain_scores);

  bool DEBUG = false;
  if (DEBUG) {
    string outfile = "shuffled_scores_tmp.txt";
    std::ofstream of;
    of.open(outfile.c_str());
    std::ostream out(of.rdbuf());
    for (size_t node_id = 0; node_id < n_nodes; ++node_id)
      for (size_t i = 0; i < domain_scores[node_id].size(); ++i)
        out << node_id << "\t" << domain_scores[node_id][i] << endl;
  }

  // combine all random HYPO domain scores from all species
  for (size_t node_id = 0; node_id < n_nodes; ++node_id)
    random_scores.insert(random_scores.end(), domain_scores[node_id].begin(),
                         domain_scores[node_id].end());
  sort(random_scores.begin(), random_scores.end());
}


static void
assign_p_values(const vector<double> &random_scores,
                const vector<vector<double> >&observed_scores,
                vector<vector<double> >&p_values) {
  const double n_randoms =
    random_scores.size() == 0 ? 1 : random_scores.size();
  for (size_t node_id = 0; node_id < observed_scores.size(); ++node_id)
    for (size_t i = 0; i < observed_scores[node_id].size(); ++i) {
      double pval = (random_scores.end() -
                     upper_bound(random_scores.begin(), random_scores.end(),
                                 observed_scores[node_id][i]))/n_randoms;
      p_values[node_id].push_back(pval);
    }
}

double
get_fdr_cutoff(const vector<vector<double> > &scores, const double fdr) {
  if (fdr <= 0)
    return numeric_limits<double>::max();
  else if (fdr > 1)
    return numeric_limits<double>::min();
  vector<double> local;
  for (size_t node_id = 0; node_id < scores.size(); ++ node_id)
    local.insert(local.end(), scores[node_id].begin(),scores[node_id].end());
  std::sort(local.begin(), local.end());
  size_t i = 0;
  for (; i < local.size() - 1 &&
         local[i+1] < fdr*static_cast<double>(i+1)/local.size(); ++i);
  return local[i];
}

static void
build_domain(const vector<MSite> &sites,
             const vector<pair<size_t, size_t> > &blocks_by_sp,
             const vector<vector<double> > &meth_probs,
             const vector<string> &node_names,
             vector<vector<GenomicRegion> > &domains) {

  for (size_t node_id = 0; node_id < node_names.size(); ++node_id) {
    for (size_t i = 0; i < blocks_by_sp.size(); ++i) {
      size_t start = blocks_by_sp[i].first;
      size_t end = blocks_by_sp[i].second;
      vector<bool> classes;
      for (size_t pos = start; pos <= end; ++pos) {
        double p = meth_probs[pos][node_id];
        bool state = p > 0.5;
        if (pos == start) {
          classes.push_back(state);
          if (!state) {
            // region score is # of sites
            domains[node_id].push_back(GenomicRegion(sites[pos].chrom,
                                                     sites[pos].pos,
                                                     sites[pos].pos + 1,
                                                     node_names[node_id],
                                                     1.0, '+'));
          }
        } else {
          if (classes.back() == state) {
            if (!state) {
              domains[node_id].back().set_end(sites[pos].pos + 1);
              const double sc = domains[node_id].back().get_score() + 1;
              domains[node_id].back().set_score(sc);
            }
          } else {
            classes.push_back(state);
            if (!state) {
              domains[node_id].push_back(GenomicRegion(sites[pos].chrom,
                                                       sites[pos].pos,
                                                       sites[pos].pos + 1,
                                                       node_names[node_id],
                                                       1.0, '+'));
            }
          }
        }
      }
    }
  }
}

// blocks for individual species
static void
build_domain(const vector<MSite> &sites,
             const vector<vector<pair<size_t, size_t> > > &blocks_by_sp,
             const vector<vector<double> > &meth_probs,
             const vector<string> &node_names,
             vector<vector<GenomicRegion> > &domains) {
  const size_t n_nodes = node_names.size();
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    for (size_t i = 0; i < blocks_by_sp[node_id].size(); ++i) {
      size_t start = blocks_by_sp[node_id][i].first;
      size_t end = blocks_by_sp[node_id][i].second;
      vector<bool> classes;
      for (size_t pos = start; pos <= end; ++pos) {
        double p = meth_probs[pos][node_id];
        bool state = p > 0.5;
        if (pos == start) {
          classes.push_back(state);
          if (!state) {
            domains[node_id].push_back(GenomicRegion(sites[pos].chrom,
                                                     sites[pos].pos,
                                                     sites[pos].pos + 1,
                                                     node_names[node_id],
                                                     1.0, '+'));
          }
        } else {
          if (classes.back() == state) {
            if (!state) {
              domains[node_id].back().set_end(sites[pos].pos + 1);
              const double sc = domains[node_id].back().get_score() + 1;
              domains[node_id].back().set_score(sc);
            }
          } else {
            classes.push_back(state);
            if (!state) {
              domains[node_id].push_back(GenomicRegion(sites[pos].chrom,
                                                       sites[pos].pos,
                                                       sites[pos].pos + 1,
                                                       node_names[node_id],
                                                       1.0, '+'));
            }
          }
        }
      }
    }
  }
}


static bool
has_missing_site(const vector<vector<double> > &tree_probs) {
  for (size_t i= 0; i < tree_probs.size(); ++i)
    for (size_t j = 0; j < tree_probs[0].size(); ++j)
      if (missing_meth_value(tree_probs[i][j]))
        return true;

  return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////          MAIN             ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    size_t rng_seed = numeric_limits<size_t>::max();

    string outfile;

    size_t desert_size = 1000;

    // run mode flags
    bool VERBOSE = false;
    double fdr_cutoff = numeric_limits<double>::max();
    double fdr = 0.01;

    /********************* COMMAND LINE OPTIONS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), "segment epigenomic states for "
                           "species with specified evolutionary relationships",
                           "<param-file> <meth-table>");
    opt_parse.add_opt("seed", 's', "rng seed (default: none)",
                      false, rng_seed);
    opt_parse.add_opt("desert", 'd', "desert size (default: " +
                      to_string(desert_size) + ")",
                      false, desert_size);
    opt_parse.add_opt("fdr", 'f', "FDR cutoff (default: " +
                      to_string(fdr) + ")",
                      false, fdr);
    opt_parse.add_opt("verbose", 'v', "print more run info "
                      "(default: " + string(VERBOSE ? "true" : "false") + ")",
                      false, VERBOSE);
    opt_parse.add_opt("out", 'o', "output file name (default: stdout)",
                      true, outfile);

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
    const string paramfile = leftover_args.front();
    const string meth_table_file = leftover_args.back();
    /******************** END COMMAND LINE OPTIONS ********************/

    /******************** INITIALIZE PARAMETERS *******************************/
    PhyloTreePreorder t;
    param_set params;
    params.read(paramfile, t);
    if (VERBOSE)
      cerr << "[given params={" << params << "}]" << endl;

    /******************** LOAD PHYLOGENETIC TREE ******************************/
    vector<size_t> subtree_sizes;
    t.get_subtree_sizes(subtree_sizes);
    vector<double> branches;
    t.get_branch_lengths(branches);

    t.assign_missing_node_names();
    vector<string> node_names;
    t.get_node_names(node_names);

    const size_t n_nodes = subtree_sizes.size();

    vector<size_t> parent_ids;
    get_parent_id(subtree_sizes, parent_ids);

    if (VERBOSE)
      cerr << "[tree:]\n" << t.tostring() << endl;

    /******************* READ THE METHYLATION DATA ****************************/
    if (VERBOSE)
      cerr << "[reading methylation data (states or posteriors)]" << endl;

    vector<MSite> sites;
    vector<vector<double> > tree_probs;
    vector<string> meth_table_species;
    read_meth_table(meth_table_file, sites, meth_table_species, tree_probs);

    if (meth_table_species.size() != n_nodes) {
      cerr << "inconsistent tree sizes between:" << endl
           << meth_table_file << endl
           << paramfile << endl;
      return EXIT_SUCCESS;
    }

    bool use_marks = has_missing_site(tree_probs);

    vector<pair<size_t, size_t> > blocks;
    vector<vector<bool> > marks;
    vector<vector<pair<size_t, size_t> > > blocks_by_sp;
    bool DEBUG = false;

    if (!use_marks) {
      if (VERBOSE)
        cerr << "[separating deserts]" << endl;
      separate_regions(desert_size, sites, blocks);
      if (VERBOSE)
        cerr << "number of blocks: " << blocks.size() << endl;

      /******* Get random domain score and set p-values ********/
      if (VERBOSE)
        cerr << "---- shuffle cpgs ---" << endl;
      vector<double> shuffled_domain_scores(n_nodes);
      shuffle_cpgs(rng_seed, subtree_sizes, blocks, tree_probs,
                   shuffled_domain_scores);

      if (VERBOSE)
        cerr << "---- compute domain scores ---" << endl;
      vector<vector<double> > domain_scores(n_nodes);
      get_domain_scores(subtree_sizes, tree_probs, blocks, domain_scores);

      if (VERBOSE)
        cerr << "---- Assign p-values ---" << endl;
      vector<vector<double> > p_values(n_nodes);
      assign_p_values(shuffled_domain_scores, domain_scores, p_values);

      if (DEBUG) {
        string tmpoutfile = "scores.txt";
        std::ofstream of;
        of.open(tmpoutfile.c_str());
        std::ostream out(of.rdbuf());
        for (size_t node_id = 0; node_id < n_nodes; ++node_id)
          for (size_t i = 0; i < domain_scores[node_id].size(); ++i)
            out << node_id << "\t" << domain_scores[node_id][i] << "\t"
                << p_values[node_id][i] << endl;
      }

      if (fdr_cutoff == numeric_limits<double>::max())
        fdr_cutoff = get_fdr_cutoff(p_values, fdr);

      if (VERBOSE)
        cerr << "[Building HYPO domains in each species]" << endl;
      vector<vector<GenomicRegion> > domains(n_nodes);
      build_domain(sites, blocks, tree_probs, node_names, domains);

      /***** output the domains ******/
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      if (!out)
        throw SMITHLABException("bad output file: " + outfile);
      for (size_t node_id = 0; node_id < domains.size(); ++node_id) {
        size_t good_hmr_count = 0;
        for (size_t i = 0; i < domains[node_id].size(); ++i)
          if (p_values[node_id][i] < fdr_cutoff) {
            domains[node_id][i].set_name(node_names[node_id] + "_HYPO"+
                                         smithlab::toa(good_hmr_count++));
            out << domains[node_id][i] << '\n';
          }
      }

    } else {

      if (VERBOSE)
        cerr << "[marking sites]" << endl;
      mark_useable_sites(subtree_sizes, sites, tree_probs, desert_size, marks);
      marks_to_blocks(desert_size, subtree_sizes, sites, marks, blocks_by_sp);
      if (VERBOSE) {
        vector<size_t> counts(n_nodes);
        for (size_t i = 0; i < marks.size(); ++i)
          for (size_t j = 0; j < n_nodes; ++j)
            if (marks[i][j]) ++counts[j];

        for (size_t j = 0; j < n_nodes; ++j)
          cerr << "node " << j << " has " << counts[j]
               << " useable sites and "
               << blocks_by_sp[j].size() << " blocks"<< endl;
      }

      /******* Get random domain score and set p-values ********/
      vector<double> shuffled_domain_scores(n_nodes);
      shuffle_cpgs(rng_seed, subtree_sizes, blocks_by_sp, tree_probs,
                   shuffled_domain_scores);

      vector<vector<double> > domain_scores(n_nodes);
      get_domain_scores(subtree_sizes, tree_probs, blocks_by_sp, domain_scores);

      vector<vector<double> > p_values(n_nodes);
      assign_p_values(shuffled_domain_scores, domain_scores, p_values);

      if (fdr_cutoff == numeric_limits<double>::max())
        fdr_cutoff = get_fdr_cutoff(p_values, fdr);

      if (DEBUG) {
        string tmpoutfile = "scores.txt";
        std::ofstream of;
        of.open(tmpoutfile.c_str());
        std::ostream out(of.rdbuf());
        for (size_t node_id = 0; node_id < n_nodes; ++node_id)
          for (size_t i = 0; i < domain_scores[node_id].size(); ++i)
            out << node_id << "\t" << domain_scores[node_id][i] << "\t"
                << p_values[node_id][i] << endl;
      }

      if (VERBOSE)
        cerr << "[Building HYPO domains in each species]" << endl;
      vector<vector<GenomicRegion> > domains(n_nodes);
      build_domain(sites, blocks_by_sp, tree_probs, node_names, domains);

      /***** output the domains ******/
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      if (!out)
        throw SMITHLABException("bad output file: " + outfile);
      for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
        size_t good_hmr_count = 0;
        for (size_t i = 0; i < domains[node_id].size(); ++i)
          if (p_values[node_id][i] < fdr_cutoff) {
            domains[node_id][i].set_name(node_names[node_id] + "_HYPO"+
                                         smithlab::toa(good_hmr_count++));
            out << domains[node_id][i] << '\n';
          }
      }

    }
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
