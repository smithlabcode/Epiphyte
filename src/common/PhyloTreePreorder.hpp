/*    Copyright (C) 2015 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith and Jenny Qu
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

#ifndef PHYLOTREE_PREORDER_HPP
#define PHYLOTREE_PREORDER_HPP

#include <string>
#include <vector>

#include "PhyloTree.hpp"

using std::string;
using std::vector;

class PhyloTreePreorder : public PhyloTree {
public:
  void
  get_subtree_sizes(vector<size_t> &subtree_sizes) const {
    get_subtree_sizes(root, subtree_sizes);
  }
  void
  get_branch_lengths(vector<double> &branch_lengths) const {
    get_branch_lengths(root, branch_lengths);
  }
  void
  get_leaf_names(vector<string> &names) const {
    get_leaf_names(root, names);
  };
  void
  get_node_names(vector<string> &names) const {
    get_node_names(root, names);
  };

  void
  assign_missing_node_names() {
    assign_missing_node_names(root, 0);
  }
  void
  set_branch_lengths(vector<double> branch_lengths) {
    set_branch_lengths(root, branch_lengths);
  }

private:
  void
  get_subtree_sizes(const PhyloTree::PTNode &node,
                    vector<size_t> &subtree_sizes) const;
  void
  get_branch_lengths(const PhyloTree::PTNode &node,
                     vector<double> &branch_lengths) const;
  void
  get_leaf_names(const PhyloTree::PTNode &node,
                 vector<string> &names) const;
  void
  get_node_names(const PhyloTree::PTNode &node,
                 vector<string> &names) const;
  void
  set_branch_lengths(PhyloTree::PTNode &node,
                     vector<double> &branch_lengths);
  void
  assign_missing_node_names(PhyloTree::PTNode &node,
                            size_t count);
};

void
get_parent_id(const std::vector<size_t> &subtree_sizes,
              std::vector<size_t> &parent_id);

bool
has_same_species_order(const PhyloTreePreorder &the_tree,
                       const std::vector<std::string> &species_names);

inline bool
is_root(const size_t node_id) {return node_id == 0;}

inline bool
is_leaf(const size_t subtree_size) {return subtree_size == 1;}

inline bool
is_binary(const vector<size_t> &subtree_sizes) {
  // ADS: this function seems not to be used
  return (subtree_sizes[0] == 1 + subtree_sizes[1] +
          subtree_sizes[subtree_sizes[1] + 1]);
}

void
subtree_sizes_to_leaves_preorder(const std::vector<size_t> &subtree_sizes,
                                 std::vector<size_t> &leaves_preorder);

void
get_degrees(const std::vector<size_t> &subtree_sizes, std::vector<size_t> &degrees);

bool
is_semi_binary(const std::vector<size_t> &degrees);

size_t
count_leaves(const std::vector<size_t> &subtree_sizes);

void
get_children(const size_t node_id, const std::vector<size_t> &subtree_sizes,
             std::vector<size_t> &children);

size_t
leafsize(const std::vector<size_t> &subtree_sizes);
#endif
