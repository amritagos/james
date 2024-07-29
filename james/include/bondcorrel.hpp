#pragma once

#include "network_base.hpp"
#include "undirected_network.hpp"
#include <cstddef>
#include <optional>
#include <vector>

namespace James::Bond::Correlation {
/*Creates a flattened vector containing a binary value for each pair of atoms
 * (i,j), corresponding to whether there is a connection [c(i,j)=1] or not
 * [c(i,j)=0]. The values are arranged in the order (0,1), (0,2)...(1,2),
 * (1,3)... such that i<j (equal values ommitted since the bond cannot be
 * between the atom itself).Takes the UndirectedNetwork as an
 * input. */
template <typename WeightType = double>
std::vector<int> bond_connection_info_at_tau(
    const Graph::UndirectedNetwork<WeightType> &network) {
  std::vector<int> c_ij{}; // flattened vector of c(i,j) pairs. 0 if there is no
                           // connection and 1 if there is a connection
  const size_t n_atoms = network.n_agents();
  int c_val{};

  // Loop through i j pairs
  for (size_t i = 0; i < n_atoms - 1; i++) {
    for (size_t j = i + 1; j < n_atoms; j++) {
      // For a particular pair (i,j)
      if (network.connection_exists(i, j)) {
        c_val = 1;
      } else {
        c_val = 0;
      }
      c_ij.push_back(c_val);
    }
  }

  return c_ij;
}
} // namespace James::Bond::Correlation