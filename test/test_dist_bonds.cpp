#include "bondfinder.hpp"
#include "pairtypes.hpp"
#include "system.hpp"
#include "undirected_network.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <vector>

TEST_CASE("Test that distance-based bonds can be created by providing pairwise "
          "cutoffs or a global cutoff",
          "[DistanceBonds]") {
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
  size_t o_index = 0;
  auto pair_oh = James::Bond::Pair(o_type, h_type);

  // The undirected network
  auto network = Graph::UndirectedNetwork<double>(system.n_atoms());

  // Create distance based cutoffs (intra-molecular water bonds)
  auto pairs = std::vector<James::Bond::Pair>{pair_oh};
  auto cutoffs = std::vector<double>{1.0};
  James::Bond::add_distance_based_bonds(network, system, pairs, cutoffs);

  // There should be two bonds between the O and H
  REQUIRE(network.n_edges(o_index) == 2);

  // The total number of edges should also be 2
  REQUIRE(network.n_edges() == 2);

  // Test the global cutoff distance criterion
  network.clear();
  James::Bond::add_distance_based_bonds(network, system, 3.5);

  // The total number of edges should be 5
  REQUIRE(network.n_edges() == 5);
}