#include "bondfinder.hpp"
#include "fmt/core.h"
#include "pairtypes.hpp"
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
}