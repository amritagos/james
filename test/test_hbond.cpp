#include "bondcorrel.hpp"
#include "bondfinder.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "fmt/core.h"
#include "fmt/ostream.h"
#include "system.hpp"
#include "undirected_network.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <vector>

TEST_CASE("Test that a hydrogen bond can be found between a Cl- ion and water "
          "molecule",
          "[Hbond]") {
  // Create a system with one water molecule and one Chloride ion
  auto system = James::Atoms::System();
  system.box = std::vector<double>{49.4752565268, 49.4752565268, 49.4752565268};
  system.boxLo =
      std::vector<double>{-12.457628299, -12.457628299, -12.457628299};
  // Positions of all the atoms (O, H, H, Cl)
  auto positions =
      std::vector<std::vector<double>>{{26.6679, 23.3139, 22.8734},
                                       {26.0598, 22.5748, 22.8577},
                                       {26.7405, 23.5454, 23.7993},
                                       {25.1689, 20.8364, 22.0004}};

  system.push_back(James::Atoms::Atom(1, 1, 1, positions[0])); // O atom
  system.push_back(James::Atoms::Atom(2, 2, 1, positions[1])); // H atom (h1)
  system.push_back(James::Atoms::Atom(3, 2, 1, positions[2])); // H atom (h2)
  system.push_back(James::Atoms::Atom(4, 4, 2, positions[3])); // Cl ion

  // Type and index information of all the atoms
  int o_type = 1;
  int h_type = 2;
  int cl_type = 4;
  size_t cl_index = 3;

  // The undirected network
  auto network = Graph::UndirectedNetwork<double>(system.n_atoms());

  // Information needed for hydrogen bonds
  auto donor_atom_types = std::vector<int>{o_type};
  auto acceptor_atom_types = std::vector<int>{cl_type, o_type};
  auto h_atom_types = std::vector<int>{h_type};
  double donor_acceptor_cutoff = 3.2;
  double max_angle_deg = 30;
  // Add the hydrogen bond between the H and the Cl
  James::Bond::add_hbonds(network, system, donor_atom_types,
                          acceptor_atom_types, h_atom_types,
                          donor_acceptor_cutoff, max_angle_deg, false);

  // There should be one hydrogen bond between the Cl and H
  REQUIRE(network.n_edges(cl_index) == 1);

  // The total number of edges should also be 1
  REQUIRE(network.n_edges() == 1);

  // Find the H in the hydrogen bond
  auto h_neighbours = network.get_neighbours(cl_index);
  INFO(fmt::format("Hydrogen index in the hydrogen bond is {} \n",
                   h_neighbours[0]));

  // Get the angle HDA and check that it is smaller than the maximum angle
  // cutoff A : acceptor (Cl-) D: donor (O) You can also get all the positions
  // using system.
  auto a_pos = system.atoms[3].position; // Cl position
  auto d_pos = system.atoms[0].position; // O position
  for (auto h_index : h_neighbours) {
    auto h_pos = system.atoms[h_index].position;
    double hda_angle =
        James::Misc::angleABCdeg(h_pos, d_pos, a_pos, system.box);
    INFO(fmt::format("HDA angle in degrees = {}\n", hda_angle));
    REQUIRE(hda_angle < 30);
  }

  // Check that if you create a flattened array containing the hydrogen bond
  // information of dissimilar (i,j) pairs, it is correct
  auto c_ij_expected = std::vector<int>{0, 0, 0, 0, 1, 0};
  auto c_ij = James::Bond::Correlation::bond_connection_info_at_tau(network);

  REQUIRE_THAT(c_ij, Catch::Matchers::RangeEquals(c_ij_expected));

  // Should work for multiple networks too
  auto networks =
      std::vector<Graph::UndirectedNetwork<double>>{network, network};
  auto c_ij_multiple_expected =
      std::vector<std::vector<int>>{c_ij_expected, c_ij_expected};
  auto c_ij_multiple =
      James::Bond::Correlation::bond_connection_info_time_series(networks,
                                                                 false, 1);

  REQUIRE_THAT(c_ij_multiple,
               Catch::Matchers::RangeEquals(c_ij_multiple_expected));

  // Multiple networks using the continuous hydrogen bond definition
  auto network_empty = Graph::UndirectedNetwork<double>(system.n_atoms());
  // Make the first network empty
  networks[0] = network_empty;
  c_ij_multiple_expected =
      std::vector<std::vector<int>>{{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
  // The flattened array will just be empty since the first network was empty
  // In the definition of the continuous hydrogen bond, bonds are considered
  // broken even if reformed later.
  c_ij_multiple = James::Bond::Correlation::bond_connection_info_time_series(
      networks, true, 1);

  REQUIRE_THAT(c_ij_multiple,
               Catch::Matchers::RangeEquals(c_ij_multiple_expected));
}