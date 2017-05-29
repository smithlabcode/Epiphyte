/*    Copyright (C) 2015 University of Southern California and
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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
#include <iomanip> // std::setw
#include <fstream>
#include <numeric> //std::accumulate
#include <unordered_map>
#include <random>
#include <algorithm>
#include <cmath>    //std::abs, floor
#include <limits>   //std::numeric_limits
#include <iterator>     // std::distance
#include <unistd.h>

/* from smithlab_cpp */
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

/* from methpipe */
#include "MethpipeFiles.hpp"

#include "BetaBin.hpp"

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

using std::string;
using std::vector;

using std::pair;
using std::make_pair;
using std::unordered_map;

using std::abs;
using std::numeric_limits;
using std::accumulate;
using std::inner_product;
using std::min;
using std::max;
using std::isfinite;

using std::endl;
using std::cerr;
using std::cout;
using std::istringstream;
using std::istream_iterator;
using std::setprecision;
using std::setw;



template <class T> T
btwn_min_max(const T &val, const T &the_min, const T &the_max) {
  const T minval = max(val, the_min);
  return min(minval, the_max);
}



struct Site {
  string chrom;
  size_t pos;
  char strand;
  string name;

  Site() {}
  Site(const string &chr, const size_t &position) :
    chrom(chr), pos(position) {}
};


std::ostream&
operator<<(std::ostream &out, const Site &s) {
  return out << s.chrom << '\t' << s.pos << '\t' << s.strand << '\t' << s.name;
}


static std::istream &
read_site(std::istream &in, vector<Site> &sites,
          vector<pair<double, double> > &counts) {
  Site s;
  double meth = 0.0, coverage = 0.0;
  if (in >> s.chrom >> s.pos >> s.strand >> s.name >> meth >> coverage) {
    const double meth_count = round(meth*coverage);
    counts.push_back(make_pair(meth_count, coverage - meth_count));
    sites.push_back(s);
  }
  return in;
}


static std::istream &
read_site(std::istream &in, Site &s, double &meth, size_t &coverage) {
  return (in >> s.chrom >> s.pos >> s.strand >> s.name >> meth >> coverage);
}


double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


static double
expectation_step(const vector<pair<double, double> > &counts,
                 const double &mixing,
                 const betabin &low_distro,
                 const betabin &high_distro,
                 vector<double> &low_probs,
                 vector<double> &high_probs) {

  double score = 0;

  const double low_log_mixing = log(mixing);
  assert(isfinite(low_log_mixing));
  const double high_log_mixing = log(1.0 - mixing);
  assert(isfinite(high_log_mixing));

  for (size_t i = 0; i < counts.size(); ++i) {

    const double low_part = low_log_mixing + low_distro(counts[i]);
    assert(isfinite(low_part));

    const double high_part = high_log_mixing + high_distro(counts[i]);
    assert(isfinite(high_part));

    const double denom = log_sum_log(low_part, high_part);
    assert(isfinite(denom));

    low_probs[i] = exp(low_part - denom);
    high_probs[i] = exp(high_part - denom);
    score += denom;

  }
  return score;
}



static void
maximization_step(const vector<double> &vals_a,
                  const vector<double> &vals_b,
                  const vector<double> &low_probs,
                  const vector<double> &high_probs,
                  double &mixing,
                  betabin &low_distro, betabin &high_distro) {

  low_distro.fit(vals_a, vals_b, low_probs);
  high_distro.fit(vals_a, vals_b, high_probs);

  vector<double> log_low_probs(low_probs);
  vector<double> log_high_probs(high_probs);
  for (size_t i = 0; i < log_low_probs.size(); ++i) {
    log_low_probs[i] = log(log_low_probs[i]);
    log_high_probs[i] = log(log_high_probs[i]);
  }

  mixing = smithlab::log_sum_log_vec(log_low_probs, log_low_probs.size());
  const double high_mixing =
    smithlab::log_sum_log_vec(log_high_probs, log_high_probs.size());
  const double mix_sum = ((mixing > high_mixing) ?
                          mixing + log(1 + exp(high_mixing - mixing)) :
                          high_mixing + log(1 + exp(mixing - high_mixing)));
  mixing = exp(mixing - mix_sum);
}



static void
two_state_em(const bool VERBOSE,
             const double tolerance, const size_t max_iter,
             const vector<pair<double, double> > &counts,
             vector<double> &posteriors,
             betabin &low_distro, betabin &high_distro,
             double &mixing) {

  // vals_a and vals_b are helper vectors for fitting the distributions
  vector<double> vals_a(counts.size()), vals_b(counts.size());
  for (size_t i = 0; i < counts.size(); ++i) {
    const double proportion =
      btwn_min_max(counts[i].first/(counts[i].first + counts[i].second),
                   1e-03, 1.0 - 1e-03);
    vals_a[i] = log(proportion);
    vals_b[i] = log(1.0 - proportion);
  }

  // now do the expectation maximization
  double prev_score = std::numeric_limits<double>::max();
  vector<double> low_probs(counts.size(), 0), high_probs(counts.size(), 0);

  if (VERBOSE)
    cerr << "itr" << '\t'
         << "MUx" << "[low_distr]" << '\t'
         << "(1-MU)x" << "[high_distr]" << '\t'
         << "delta" << endl;

  for (size_t itr = 0; itr < max_iter; ++itr) {
    const double score = expectation_step(counts, mixing, low_distro,
                                          high_distro, low_probs, high_probs);
    maximization_step(vals_a, vals_b, low_probs, high_probs, mixing,
                      low_distro, high_distro);

    if (VERBOSE) {
      cerr << itr << '\t'
           << setprecision(2) << std::fixed
           << mixing << ' '
           << '[' << low_distro.tostring() << ']' << '\t'
           << 1.0 - mixing << ' '
           << '[' << high_distro.tostring() << ']' << '\t'
           << setprecision(1) << std::scientific
           << (prev_score - score)/prev_score
           << endl;
    }
    if ((prev_score - score)/prev_score < tolerance)
      break;
    prev_score = score;
  }

  posteriors.swap(high_probs);
}



static void
get_posteriors(const vector<pair<double, double> > &counts,
               const double &mixing,
               const betabin &low_distro,
               const betabin &high_distro,
               vector<double> &posteriors) {
  vector<double> low_probs(counts.size(), 0), high_probs(counts.size(), 0);
  expectation_step(counts, mixing,
                   low_distro, high_distro, low_probs, high_probs);
  posteriors.swap(high_probs);
}



static void
read_parameters(const string &infile, double &mixing,
                double &low_alpha, double &low_beta,
                double &high_alpha, double &high_beta) {
  std::ifstream in(infile.c_str());
  if (!in)
    throw SMITHLABException("could not open file: " + infile);

  string label;
  if (!(in >> label >> mixing
        >> label >> low_alpha
        >> label >> low_beta
        >> label >> high_alpha
        >> label >> high_beta))
    throw SMITHLABException("bad format for params file: " + infile);
}



static void
write_parameters(const double mixing,
                 const betabin &low_distro, const betabin &high_distro,
                 const string &outfile) {

  std::ofstream out(outfile.c_str());
  if (!out)
    cerr << "bad parameter outfile: " << outfile << endl;
  else {
    out.precision(30);
    out << "MIXING\t" << mixing << endl
        << "LOW_ALPHA\t" << low_distro.alpha << endl
        << "LOW_BETA\t" << low_distro.beta << endl
        << "HIGH_ALPHA\t" << high_distro.alpha << endl
        << "HIGH_BETA\t" << high_distro.beta << endl;
  }
}



int
main(int argc, const char **argv) {

  try {

    // run mode flags
    bool VERBOSE = false;

    double TOLERANCE = 1e-10;
    size_t MAXITER = 50;

    string outfile;
    string params_in_file;
    string params_out_file;

    /************************* COMMAND LINE OPTIONS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "posteriors for "
                           "hyper-methylation state assuming independent sites",
                           "<methcounts>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("params-in", '\0', "parameters file (no training)",
                      false, params_in_file);
    opt_parse.add_opt("params-out", '\0', "write fit parameters to this file",
                      false, params_out_file);
    opt_parse.add_opt("verbose", 'v', "print more run info (default: false)",
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (outfile.empty() && params_out_file.empty()) {
      cerr << "provide at least one of "
           << "params output file or posteriors output file" << endl;
      return EXIT_SUCCESS;
    }
    if (!params_out_file.empty() && !params_in_file.empty()) {
      cerr << "provide at most one of "
           << "params input file and posteriors output file" << endl;
      return EXIT_SUCCESS;
    }
    const string meth_table_file = leftover_args.front();
    /**************************************************************************/

    std::ifstream in(meth_table_file.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + meth_table_file);

    if (VERBOSE)
      cerr << "loading data" << endl;

    vector<pair<double, double> > counts;
    vector<Site> sites;
    while (read_site(in, sites, counts)) {}
    in.close();
    if (VERBOSE)
      cerr << "total sites read: " << counts.size() << endl;

    // remove not-covered sites
    size_t j = 0;
    for (size_t i = 0; i < counts.size(); ++i) {
      if (counts[i].first + counts[i].second > 0) {
        std::swap(counts[j], counts[i]);
        std::swap(sites[j], sites[i]);
        ++j;
      }
    }
    counts.erase(counts.begin() + j, counts.end());
    sites.erase(sites.begin() + j, sites.end());
    if (VERBOSE)
      cerr << "sites with data: " << counts.size() << endl;

    // add pseudocounts
    for (size_t i = 0; i < counts.size(); ++i) {
      counts[i].first++;
      counts[i].second++;
    }

    // declare the parameters that will either be trained or specified
    // in a file
    double low_alpha = 0.0, low_beta = 0.0;
    double high_alpha = 0.0, high_beta = 0.0;
    double mixing = 0.0;

    // read the parameters if specified
    if (!params_in_file.empty()) {
      read_parameters(params_in_file, mixing,
                      low_alpha, low_beta, high_alpha, high_beta);
    }
    else {
      // initialize the Beta Binomial parameters (will be trained)
      low_alpha = 0.5;
      low_beta = 4; // low mean
      high_alpha = 4;
      high_beta = 0.5; // high mean
      mixing = 0.5; // initialize the mixing proportion to fifty-fifty
    }

    betabin low_distro(low_alpha, low_beta);
    betabin high_distro(high_alpha, high_beta);

    vector<double> posteriors;

    if (params_in_file.empty()) {
      if (VERBOSE)
        cerr << "training parameters" << endl;
      two_state_em(VERBOSE, TOLERANCE, MAXITER, counts, posteriors,
                   low_distro, high_distro, mixing);
      // write parameters if requested
      if (!params_out_file.empty())
        write_parameters(mixing, low_distro, high_distro, params_out_file);
    }
    else {
      if (VERBOSE)
        cerr << "computing posteriors" << endl;
      get_posteriors(counts, mixing, low_distro, high_distro, posteriors);
    }


    if (!outfile.empty()) {
      std::ofstream out(outfile.c_str());
      if (!out)
        throw SMITHLABException("bad output file: " + outfile);

      in.clear();
      in.open(meth_table_file.c_str());
      Site s;
      double meth = 0.0;
      size_t n_reads = 0;
      size_t j = 0;
      while (read_site(in, s, meth, n_reads)) {
        out << s << '\t' << ((n_reads > 0) ? posteriors[j] : 0.0) << '\t'
            << n_reads << endl;
        j += (n_reads > 0);
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
