#pragma once
#include "network_base.hpp"
#include "pairtypes.hpp"
#include "system.hpp"
#include "util.hpp"
#include <cstddef>
#include <vector>

namespace James::Bond {

/*
Adds hydrogen bonds, according to a geometric criterion. A donor (D) is an
electronegative atom, which is bonded to at least one hydrogen. An acceptor (A)
is an electronegative atom which accepts the hydrogen (H). The criterion is as
follows: 1) the distance DA between the donor and acceptor should be less than
or equal to 3.5 Angstroms, which corresponds to the first minimum in the O-O RDF
of water, and 2) the angle HDA should be less than or equal to 30 degrees.

If you set ignore_hydrogens to true (default), then a bond is added between the
donor and acceptor only If this is set to false, then the bond between the
hydrogen and the acceptor is added.

Adds hydrogen bonds to a given UndirectedNetwork or DirectedNetwork. For
molecular systems, you should use UndirectedNetwork.
*/
template <typename WeightType = double>
void add_hbonds(Graph::NetworkBase<WeightType> &network,
                const James::Atoms::System &system,
                const std::vector<int> &donor_atom_types,
                const std::vector<int> &acceptor_atom_types,
                const std::vector<int> &h_atom_types,
                double cutoff_distance = 3.2, double max_angle_deg = 30,
                bool ignore_hydrogens = true) {
  // Lambda for checking if a particular atom type (target) is one of allowed
  // atom types
  auto atom_type_found = [&](const std::vector<int> &atom_types, int target) {
    return std::find(atom_types.begin(), atom_types.end(), target) !=
           atom_types.end();
  };
  // Loop through all atoms
  for (size_t i_a = 0; i_a < system.atoms.size(); i_a++) {
    // Only do stuff if an acceptor has been found
    if (atom_type_found(acceptor_atom_types, system.atoms[i_a].type)) {
      // Now search for donors
      for (size_t i_d = 0; i_d < system.atoms.size(); i_d++) {
        // Skip if you're on the same atom
        if (i_a == i_d) {
          continue;
        }

        if (atom_type_found(donor_atom_types, system.atoms[i_d].type)) {
          // 1. The D-A distance should be less than or equal to 3.5 Angstrom
          // Skip if D-A distance is greater than the cutoff (3.5 Angstroms)
          if (system.distance(i_a, i_d) > cutoff_distance) {
            continue;
          }

          auto acceptor = system.atoms[i_a];
          auto donor = system.atoms[i_d];
          // Find H atoms bonded to the donor D
          // They should have the same molecule ID
          auto atom_indices_molecule =
              system.find_atoms_in_molecule(donor.mol_id.value());

          // Loop through these atoms and find H atoms
          for (auto &index : atom_indices_molecule) {
            if (atom_type_found(h_atom_types, system.atoms[index].type)) {
              // H atom found
              auto h_atom = system.atoms[index];

              // 2. The HDA angle should be less than or equal to 30 degrees
              double hda_angle =
                  Misc::angleABCdeg(h_atom.position, donor.position,
                                    acceptor.position, system.box);

              if (hda_angle > max_angle_deg) {
                continue;
              }

              // Hydrogen bond has been found between the donor and acceptor
              // Add the hydrogen bond
              if (ignore_hydrogens == true) {
                // Add a hydrogen bond between the donor and acceptor
                network.push_back_neighbour_and_weight(i_d, i_a, 1.0);
              } else {
                // Add a bond between H and the acceptor
                network.push_back_neighbour_and_weight(index, i_a, 1.0);
              }
            }
          }
        }
      }
    }
  }
}

/*
Adds bonds based on a distance based criterion.
An unordered_map of pair types and cutoffs will be created internally from the
vector of Pair objects and corresponding cutoffs. For any pair of types not
provided, no bond will be created.
*/
template <typename WeightType = double>
void add_distance_based_bonds(Graph::NetworkBase<WeightType> &network,
                              const James::Atoms::System &system,
                              std::vector<Pair> &pairs,
                              std::vector<double> &cutoffs) {
  // The number of Pair objects should correspond to an equal number of cutoff
  // values
  if (pairs.size() != cutoffs.size()) {
    std::runtime_error("Inconsistent Pair and cutoff values\n");
  }
  auto pair_cutoff_map = create_pairtype_cutoffs(
      pairs, cutoffs); // create an unordered_map with Pair objects as keys and
                       // cutoffs as the values

  // Loop through all i,j pairs of atoms
  for (size_t i_idx = 0; i_idx < system.atoms.size() - 1; i_idx++) {
    for (size_t j_idx = i_idx + 1; j_idx < system.atoms.size(); j_idx++) {
      int i_type = system.atoms[i_idx].type;
      int j_type = system.atoms[j_idx].type;
      // If the key with the Pair i_type, j_type does not exist, then skip
      if (pair_cutoff_map.contains(Pair(i_type, j_type == false))) {
        continue;
      }
      // If the Pair exists, then the cutoff can be accessed
      double cutoff = pair_cutoff_map[Pair(i_type, j_type)];
      // Get the distance
      double r_ij = system.distance(i_idx, j_idx);

      // Check that the distance between the atoms is within the cutoff
      if (r_ij <= cutoff) {
        // Add the bond
        network.push_back_neighbour_and_weight(i_idx, j_idx, 1.0);
      }
    }
  }
}

/*
Adds bonds based on a distance based criterion, using a global cutoff value and
ignoring pair types.
*/
template <typename WeightType = double>
void add_distance_based_bonds(Graph::NetworkBase<WeightType> &network,
                              const James::Atoms::System &system,
                              double global_cutoff = 3.5) {

  // Loop through all i,j pairs of atoms
  for (size_t i_idx = 0; i_idx < system.atoms.size() - 1; i_idx++) {
    for (size_t j_idx = i_idx + 1; j_idx < system.atoms.size(); j_idx++) {
      // Get the distance
      double r_ij = system.distance(i_idx, j_idx);

      // Check that the distance between the atoms is within the cutoff
      if (r_ij <= global_cutoff) {
        // Add the bond
        network.push_back_neighbour_and_weight(i_idx, j_idx, 1.0);
      }
    }
  }
}

} // namespace James::Bond