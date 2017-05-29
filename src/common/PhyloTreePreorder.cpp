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

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "PhyloTreePreorder.hpp"

using std::string;
using std::vector;


void
PhyloTreePreorder::get_subtree_sizes(const PhyloTree::PTNode &node,
                                     vector<size_t> &subtree_sizes) const {
  subtree_sizes.push_back(1);

  if (node.has_children()) {
    for (size_t i = 0; i < node.child.size(); ++i) {
      vector<size_t> child_subtree_sizes;
      get_subtree_sizes(node.child[i], child_subtree_sizes);
      subtree_sizes.front() += child_subtree_sizes.size();
      copy(child_subtree_sizes.begin(), child_subtree_sizes.end(),
           back_inserter(subtree_sizes));
    }
  }
}


void
PhyloTreePreorder::get_branch_lengths(const PhyloTree::PTNode &node,
                                      vector<double> &branch_lengths) const {
  branch_lengths.push_back(node.branch_length);
  for (size_t i = 0; i < node.child.size(); ++i)
    get_branch_lengths(node.child[i], branch_lengths);
}


void
PhyloTreePreorder::get_leaf_names(const PhyloTree::PTNode &node,
                                  vector<string> &leafnames) const {
  if (node.child.size()==0) {
    leafnames.push_back(node.name);
  } else {
    for (size_t i = 0; i < node.child.size(); ++i)
      get_leaf_names(node.child[i], leafnames);
  }
}

void
PhyloTreePreorder::get_node_names(const PhyloTree::PTNode &node,
                                  vector<string> &node_names) const {
  node_names.push_back(node.name);
  for (size_t i = 0; i < node.child.size(); ++i)
    get_node_names(node.child[i], node_names);
}

void
PhyloTreePreorder::assign_missing_node_names(PhyloTree::PTNode &node,
                                             size_t count) {
  if (node.name.empty())
    node.name = "node_" + std::to_string(count++);
  for (size_t i = 0; i < node.child.size(); ++i)
    assign_missing_node_names(node.child[i], count);
}

void
PhyloTreePreorder::set_branch_lengths(PhyloTree::PTNode &node,
                                      vector<double> &branch_lengths) {
  node.branch_length = branch_lengths[0];
  branch_lengths.erase(branch_lengths.begin());
  for (size_t i = 0; i < node.child.size(); ++i)
    set_branch_lengths(node.child[i], branch_lengths);
}


void
get_parent_id(const vector<size_t> &subtree_sizes, vector<size_t> &parent_id) {
  parent_id = vector<size_t>(subtree_sizes.size(), 0);
  for (size_t i = 0; i < subtree_sizes.size(); ++i)
    for (size_t count = 1; count < subtree_sizes[i];) {
      const size_t child_id = i + count;
      parent_id[child_id] = i;
      count += subtree_sizes[child_id];
    }
}

bool
has_same_species_order(const PhyloTreePreorder &the_tree,
                       const vector<string> &meth_table_species) {
  vector<string> leaf_names;
  the_tree.get_leaf_names(leaf_names);
  return leaf_names == meth_table_species;
}


void
subtree_sizes_to_leaves_preorder(const vector<size_t> &subtree_sizes,
                                 vector<size_t> &leaves_preorder) {
  for (size_t i = 0; i < subtree_sizes.size(); ++i)
    if (subtree_sizes[i] == 1)
      leaves_preorder.push_back(i);
}


static void
get_node_degrees(const vector<size_t> &subtree_sizes,
                 const size_t tree_start,
                 vector<size_t> &degrees) {

  assert(degrees.size() == subtree_sizes.size());

  if (subtree_sizes[tree_start] == 1) {
    degrees[tree_start] = 1;
  } else {
    if (tree_start > 0)
      degrees[tree_start] = 1;

    for (size_t count = 1; count < subtree_sizes[tree_start];) {
      const size_t next_child = tree_start + count;
      get_node_degrees(subtree_sizes, next_child, degrees);
      degrees[tree_start] += 1;
      count += subtree_sizes[next_child];
    }
  }
}


void
get_degrees(const vector<size_t> &subtree_sizes, vector<size_t> &degrees) {
  degrees = vector<size_t>(subtree_sizes.size(), 0);
  get_node_degrees(subtree_sizes, 0, degrees);
}


bool
is_semi_binary(const vector<size_t> &degrees) {

  if (degrees[0] < 2 || degrees[0] > 3)
    return false;

  for (size_t i = 1; i < degrees.size(); ++i )
    if (degrees[i] != 1 && degrees[i] != 3)
      return false;

  return true;
}

size_t
count_leaves(const vector<size_t> &subtree_sizes) {
  size_t n_leaf = 0;
  for (size_t i = 0; i < subtree_sizes.size(); ++i)
    n_leaf += is_leaf(subtree_sizes[i]);
  return n_leaf;
}

void
get_children(const size_t node_id, const vector<size_t> &subtree_sizes,
             vector<size_t> &children) {
  for (size_t c = 1; c < subtree_sizes[node_id]; c += subtree_sizes[node_id + c])
    children.push_back(node_id + c);
}

size_t
leafsize(const vector<size_t> &subtree_sizes) {
  size_t n_leaf = 0;
  for (size_t i = 0; i < subtree_sizes.size(); ++i) {
    if (subtree_sizes[i] == 1)
      ++n_leaf;
  }
  return n_leaf;
}
