#include "catch2/matchers/catch_matchers.hpp"
#include "system.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <vector>

TEST_CASE("Test the system", "[System]") {
  // Create a system with two atoms
  auto system = James::Atoms::System();
  system.box = std::vector<double>{49.35, 49.35, 49.35};
  system.boxLo =
      std::vector<double>{-12.457628299, -12.457628299, -12.457628299};
  auto positions = std::vector<std::vector<double>>{
      {26.6679, 23.3139, 22.8734}, {25.1689, 20.8364, 22.0004}};

  system.push_back(James::Atoms::Atom(1, 1, 1, positions[0]));
  system.push_back(James::Atoms::Atom(4, 4, 2, positions[1]));

  // Test that the distance between the O and Cl is 3.02442
  double dist_required = 3.02442;

  REQUIRE_THAT(system.distance(0, 1),
               Catch::Matchers::WithinRel(dist_required, 1e-5));

  // Add two H atoms
  auto h1_atom_pos = std::vector<double>{26.0598, 22.5748, 22.8577};
  auto h2_atom_pos = std::vector<double>{26.7405, 23.5454, 23.7993};
  system.push_back(James::Atoms::Atom(2, 2, 1, h1_atom_pos));
  system.push_back(James::Atoms::Atom(3, 2, 1, h2_atom_pos));

  // Find all indices in atoms such that the mol_id=1
  int mol_id = system.atoms[0].mol_id.value();
  auto indices_required = std::vector<size_t>{0, 2, 3};
  auto atom_indices = system.find_atoms_in_molecule(mol_id);

  REQUIRE_THAT(atom_indices, Catch::Matchers::RangeEquals(indices_required));

  // Should return an empty vector if there are no matches
  auto empty_indices = system.find_atoms_in_molecule(4);
  REQUIRE(empty_indices.size() == 0);
}