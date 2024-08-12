#include "bondfinder.hpp"
#include "catch2/catch_message.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "fmt/core.h"
#include "fmt/ostream.h"
#include "pairtypes.hpp"
#include "pathfinder.hpp"
#include "system.hpp"
#include "undirected_network.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <vector>

TEST_CASE("Test that ion pairs can be found for a system with an Fe3+ center, "
          "water molecules and chloride ions",
          "[IonPairs]") {
  // Create a system with six water molecules, three Cl- ions and one Fe3+ ion
  auto system = James::Atoms::System();
  system.box = std::vector<double>{49.4752565268, 49.4752565268, 49.4752565268};
  system.boxLo =
      std::vector<double>{-12.457628299, -12.457628299, -12.457628299};
  // Type and index information of all the atoms
  int o_type = 1;
  int h_type = 2;
  int fe_type = 3;
  int cl_type = 4;
  // Positions of all the atoms (O, H, H, Fe3+, Cl-)
  auto fe_pos = std::vector<double>{27.7128, 23.7272, 21.2095};
  auto cl_positions =
      std::vector<std::vector<double>>{{31.2487, 25.7388, 20.4993},
                                       {26.4904, 26.8173, 18.6842},
                                       {25.1689, 20.8364, 22.0004}};
  auto wat_positions = std::vector<std::vector<double>>{
      {27.7291, 25.7896, 21.3085}, {27.4666, 26.3876, 20.6086},
      {27.9001, 26.3588, 22.0589}, {26.1004, 23.7862, 19.9861},
      {25.9507, 23.1961, 19.2474}, {25.7969, 24.6397, 19.6768},
      {28.8903, 24.0058, 19.5787}, {29.634, 24.6075, 19.6117},
      {28.8589, 23.7127, 18.668},  {26.6679, 23.3139, 22.8734},
      {26.0598, 22.5748, 22.8577}, {26.7405, 23.5454, 23.7993},
      {29.4036, 23.642, 22.3353},  {29.403, 23.02, 23.0629},
      {30.0916, 24.2696, 22.5568}, {27.9556, 21.7661, 20.944},
      {27.472, 20.9433, 20.8713},  {28.8714, 21.5182, 20.8169}};
  // Fe and Cl ions
  system.push_back(James::Atoms::Atom(1, fe_type, 1, fe_pos));          // Fe3+
  system.push_back(James::Atoms::Atom(2, cl_type, 2, cl_positions[0])); // Cl-
  system.push_back(James::Atoms::Atom(3, cl_type, 3, cl_positions[1])); // Cl-
  system.push_back(James::Atoms::Atom(4, cl_type, 4, cl_positions[2])); // Cl-

  // Water molecules
  for (size_t id = 5, mol_id = 5; id <= 22; id = id + 3, mol_id++) {
    system.push_back(
        James::Atoms::Atom(id, o_type, mol_id, wat_positions[id - 5])); // O
    system.push_back(
        James::Atoms::Atom(id + 1, h_type, mol_id, wat_positions[id - 4])); // H
    system.push_back(
        James::Atoms::Atom(id + 2, h_type, mol_id, wat_positions[id - 3])); // H
  }

  // The system should have 22 atoms
  REQUIRE(system.n_atoms() == 22);

  // The undirected network
  auto network = Graph::UndirectedNetwork<double>(system.n_atoms());

  // Create the bonds from the Fe center to the water molecules (to the O atoms)
  auto pair_feo = James::Bond::Pair(fe_type, o_type);
  auto pairs = std::vector<James::Bond::Pair>{pair_feo};
  auto cutoffs = std::vector<double>{2.6};
  James::Bond::add_distance_based_bonds(network, system, pairs, cutoffs);

  // 6 bonds should have been created emanating from the Fe3+ center
  // We know the Fe ion has an index of 0
  REQUIRE(network.n_edges(0) == 6);
  // There should also only be a total of 6 bonds
  REQUIRE(network.n_edges() == 6);

  // ----------------------------------------------------------------------------
  // Don't ignore hydrogens

  // In this case, intramolecular water bonds will be needed to ensure
  // connectivity
  auto pair_oh = James::Bond::Pair(o_type, h_type);
  pairs = std::vector<James::Bond::Pair>{pair_oh};
  cutoffs = std::vector<double>{1.0};
  James::Bond::add_distance_based_bonds(network, system, pairs, cutoffs);
  // 6*2 new bonds were created
  REQUIRE(network.n_edges() == 18);

  // Information needed for hydrogen bonds
  auto donor_atom_types = std::vector<int>{o_type};
  auto acceptor_atom_types = std::vector<int>{cl_type, o_type};
  auto h_atom_types = std::vector<int>{h_type};
  double donor_acceptor_cutoff = 3.2;
  double max_angle_deg = 30;
  // Add the hydrogen bond between the H and the acceptor atom types (Cl and O)
  James::Bond::add_hbonds(network, system, donor_atom_types,
                          acceptor_atom_types, h_atom_types,
                          donor_acceptor_cutoff, max_angle_deg, false);

  // Stuff needed for finding ion pairs (Fe and Cl at the end points)
  size_t fe_index = 0; // In this case we know that Fe has an index of 0. I
                       // would loop through the atoms ordinarily I guess
  auto destination_atom_types = std::vector<int>{cl_type};
  auto intermediate_atom_types = std::vector<int>{h_type, o_type};
  int max_depth = 5;
  // Expected ion pair paths
  auto ion_pairs_with_h_required = std::vector<std::vector<int>>{
      {1, 11, 10, 0}, {2, 5, 4, 0}, {3, 14, 13, 0}};
  auto ion_pairs_with_hydrogens = James::Path::find_ion_pairs(
      fe_index, network, system, destination_atom_types,
      intermediate_atom_types, max_depth, James::Path::WriteIdentifier::Index);
  INFO(fmt::format("Number of ion pairs found is {}",
                   ion_pairs_with_hydrogens.size()));

  REQUIRE_THAT(ion_pairs_with_hydrogens,
               Catch::Matchers::RangeEquals(ion_pairs_with_h_required));
  // ----------------------------------------------------------------------------
  // Ignore hydrogens: this is probably what we will end up doing for finding
  // ion pairs
  network.clear();
  // Create the bonds from the Fe center to the water molecules (to the O atoms)
  pairs = std::vector<James::Bond::Pair>{pair_feo};
  cutoffs = std::vector<double>{2.6};
  James::Bond::add_distance_based_bonds(network, system, pairs, cutoffs);

  // No need for intramolecular water bonds
  // Add the hydrogen bond between the donors (O) and the acceptor atom types
  // (Cl and O)
  James::Bond::add_hbonds(network, system, donor_atom_types,
                          acceptor_atom_types, h_atom_types,
                          donor_acceptor_cutoff, max_angle_deg, true);

  // Find the ion pairs and check them
  auto ion_pairs_no_h_required =
      std::vector<std::vector<int>>{{1, 10, 0}, {2, 4, 0}, {3, 13, 0}};
  auto ion_pairs_no_hydrogens = James::Path::find_ion_pairs(
      fe_index, network, system, destination_atom_types,
      intermediate_atom_types, max_depth, James::Path::WriteIdentifier::Index);
  INFO(fmt::format("Number of ion pairs found (ignoring H) is {}",
                   ion_pairs_with_hydrogens.size()));

  REQUIRE_THAT(ion_pairs_no_hydrogens,
               Catch::Matchers::RangeEquals(ion_pairs_no_h_required));

  // Check that if you set max_depth to 1 (maximum length of ion pair including
  // ends would be 2), you will find zero ion pairs.
  max_depth = 1;
  ion_pairs_no_hydrogens = James::Path::find_ion_pairs(
      fe_index, network, system, destination_atom_types,
      intermediate_atom_types, max_depth);
  REQUIRE(ion_pairs_no_hydrogens.size() == 0);
}