#pragma once
#include "network_base.hpp"
#include "network_operations.hpp"
#include "system.hpp"
#include <cstddef>
#include <cstdio>
#include <vector>

namespace James::Path {

enum class WriteIdentifier { AtomID, Index };

/* Converts a path with indices in it to atom IDs from the System object */
inline void convert_path_to_ids(std::vector<int> &path,
                                const James::Atoms::System &system) {
  for (auto &ele : path) {
    ele = system.atoms[ele].id;
  }
}

/*Checks to make sure that intermediate members are only of the allowed atom
 * types. This is only an issue if the size is greater than 2, since the end
 * points are vetted (presumably) */
inline bool
check_ion_pair_path(const std::vector<int> &test_path,
                    const James::Atoms::System &system,
                    const std::vector<int> &intermediate_atom_types) {
  // Lambda for checking if a particular atom type (target) is one of allowed
  // atom types
  auto atom_type_found = [&](const std::vector<int> &atom_types, int target) {
    return std::find(atom_types.begin(), atom_types.end(), target) !=
           atom_types.end();
  };

  if (test_path.size() > 2) {
    // Check the intermediate atom types (no need to check the end points)
    for (size_t i = 1; i < test_path.size() - 1; i++) {
      size_t i_atom = test_path[i];
      if (atom_type_found(intermediate_atom_types, system.atoms[i_atom].type) ==
          false) {
        return false;
      }
    }
    return true;
  } else if (test_path.size() == 2) {
    return true;
  } else {
    return false;
  }
}

/*
Finds shortest paths, given an UndirectedNetwork or DirectedNetwork, a System
object. For molecular systems, you should use UndirectedNetwork. All shortest
paths from a given source and a given destination atom type are found.
Intermediate members can only be allowed atom types
*/
template <typename WeightType = double>
std::vector<std::vector<int>>
find_ion_pairs(size_t source, Graph::NetworkBase<WeightType> &network,
               const James::Atoms::System &system,
               const std::vector<int> &destination_atom_types,
               const std::vector<int> &intermediate_atom_types,
               std::optional<int> max_depth,
               WriteIdentifier identifier = WriteIdentifier::AtomID) {
  // Lambda for checking if a particular atom type (target) is one of allowed
  // atom types
  auto atom_type_found = [&](const std::vector<int> &atom_types, int target) {
    return std::find(atom_types.begin(), atom_types.end(), target) !=
           atom_types.end();
  };
  const auto n_atoms =
      network.n_agents(); // Number of atoms (agents) in the network and system
  // Parent list for each node in the network, for a given source
  auto parent = std::vector<std::vector<int>>{n_atoms, std::vector<int>{}};
  // Depth of each node from the given source, initialized to INT_MAX first
  auto depth_level(std::vector<int>(network.n_agents(), INT_MAX));
  // Vector of vectors containing the ion pairs
  std::vector<std::vector<int>> ion_pairs{};

  // Perform the BFS over the graph upto max_depth
  Graph::bfs(network, parent, depth_level, source, max_depth);

  // Go through atoms to find the destination atom types
  for (size_t dest_idx = 0; dest_idx < n_atoms; dest_idx++) {
    // Skip the source itself and all atoms which are not the destination type
    if (dest_idx == source ||
        atom_type_found(destination_atom_types, system.atoms[dest_idx].type) ==
            false) {
      continue;
    }
    // A destination has been found
    // If at least one path exists between the source and destination
    if (Graph::path_exists_to_destination(depth_level, dest_idx)) {
      // Obtain the paths by iterating up the parents, since we know that they
      // exist Will contain the resultant paths
      auto shortest_paths = std::vector<std::vector<int>>{};
      auto path =
          std::vector<int>{}; // Needed for recursive reconstruct_paths function
      Graph::reconstruct_paths(parent, shortest_paths, path, dest_idx);
      // Sanity check for the paths found
      for (auto &test_path : shortest_paths) {
        if (check_ion_pair_path(test_path, system, intermediate_atom_types)) {
          if (identifier == WriteIdentifier::AtomID) {
            convert_path_to_ids(test_path, system);
            ion_pairs.push_back(test_path);
          } else if (identifier == WriteIdentifier::Index) {
            ion_pairs.push_back(test_path);
          }
        }
      }
    }
  }

  return ion_pairs;
}

} // namespace James::Path