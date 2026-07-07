#include "bondfinder.hpp"
#include "pairtypes.hpp"
#include "pathfinder.hpp"
#include "system.hpp"
#include "undirected_network.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <algorithm>
#include <optional>
#include <vector>

TEST_CASE(
    "Inter-chain ion-mediated paths between charged residues can be found",
    "[ProteinIonPaths]") {
  // Atom types for the implicit biomolecular model
  const int positive_residue_type = 1;
  const int negative_residue_type = 2;
  const int na_type = 5;
  const int cl_type = 6;

  // Construct a tiny artificial system.
  //
  // Valid inter-chain path:
  //
  //   index: 0        4   5   2
  //   type:  +res -- Cl--Na-- -res
  //   mol:   1                 2
  //
  // Intra-chain path that should be filtered out:
  //
  //   index: 0        6   7   1
  //   type:  +res -- Cl--Na-- -res
  //   mol:   1                 1
  //
  // Unphysical, but geometrically possible path:
  //
  //   index: 0        8   9   10
  //   type:  +res -- Na--Cl-- -res
  //   mol:   1                  3

  auto ids =
      std::vector<int>{101, 102, 201, 202, 501, 502, 503, 504, 601, 602, 301};

  auto types = std::vector<int>{
      positive_residue_type, // 0
      negative_residue_type, // 1
      negative_residue_type, // 2
      positive_residue_type, // 3
      cl_type,               // 4
      na_type,               // 5
      cl_type,               // 6
      na_type,               // 7
      na_type,               // 8
      cl_type,               // 9
      negative_residue_type  // 10
  };

  auto mol_ids = std::optional<std::vector<int>>{
      std::vector<int>{1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 3}};

  auto positions = std::vector<std::vector<double>>{
      {0.0, 0.0, 0.0}, // 0: positive residue, chain 1

      {0.0, 3.0, 0.0}, // 1: negative residue, chain 1
      {3.0, 0.0, 0.0}, // 2: negative residue, chain 2
      {9.0, 9.0, 9.0}, // 3: positive residue, chain 2, far away

      {1.0, 0.0, 0.0}, // 4: Cl
      {2.0, 0.0, 0.0}, // 5: Na

      {0.0, 1.0, 0.0}, // 6: Cl
      {0.0, 2.0, 0.0}, // 7: Na

      {-1.0, 0.0, 0.0}, // 8: Na, weird contact with positive residue
      {-2.0, 0.0, 0.0}, // 9: Cl
      {-3.0, 0.0, 0.0}, // 10: negative residue, chain 3
  };

  auto box =
      std::optional<std::vector<double>>{std::vector<double>{20.0, 20.0, 20.0}};
  auto box_lo =
      std::optional<std::vector<double>>{std::vector<double>{0.0, 0.0, 0.0}};

  auto system =
      James::Atoms::System(ids, types, positions, mol_ids, box, box_lo);

  auto network = Graph::UndirectedNetwork<double>(system.n_atoms());

  // Allowed graph contacts.
  //
  // The first three are the physically expected contacts:
  //
  //   pos -- Cl
  //   neg -- Na
  //   Na  -- Cl
  //
  // The last two are deliberately included so that the weird path
  //
  //   pos -- Na -- Cl -- neg
  //
  // can also be detected by the pathfinder.
  auto pairs = std::vector<James::Bond::Pair>{
      James::Bond::Pair(positive_residue_type, cl_type),
      James::Bond::Pair(negative_residue_type, na_type),
      James::Bond::Pair(na_type, cl_type),

      James::Bond::Pair(positive_residue_type, na_type),
      James::Bond::Pair(negative_residue_type, cl_type),
  };

  auto cutoffs = std::vector<double>{1.1, 1.1, 1.1, 1.1, 1.1};

  James::Bond::add_distance_based_bonds(network, system, pairs, cutoffs);

  // Intended edges:
  //
  //   0 -- 4 -- 5 -- 2
  //   0 -- 6 -- 7 -- 1
  //   0 -- 8 -- 9 -- 10
  REQUIRE(network.n_edges() == 9);

  const size_t source = 0;
  auto destination_atom_types = std::vector<int>{negative_residue_type};
  auto intermediate_atom_types = std::vector<int>{na_type, cl_type};

  auto all_paths = James::Path::find_ion_pairs(
      source, network, system, destination_atom_types, intermediate_atom_types,
      std::optional<int>{4}, James::Path::WriteIdentifier::Index);

  auto expected_all_paths = std::vector<std::vector<int>>{
      {1, 7, 6, 0},
      {2, 5, 4, 0},
      {10, 9, 8, 0},
  };

  std::sort(all_paths.begin(), all_paths.end());
  std::sort(expected_all_paths.begin(), expected_all_paths.end());

  REQUIRE_THAT(all_paths, Catch::Matchers::RangeEquals(expected_all_paths));

  auto interchain_paths = James::Path::find_interchain_ion_pairs(
      source, network, system, destination_atom_types, intermediate_atom_types,
      std::optional<int>{4}, true, James::Path::WriteIdentifier::Index);

  auto expected_interchain_paths = std::vector<std::vector<int>>{
      {2, 5, 4, 0},
      {10, 9, 8, 0},
  };

  std::sort(interchain_paths.begin(), interchain_paths.end());
  std::sort(expected_interchain_paths.begin(), expected_interchain_paths.end());

  REQUIRE_THAT(interchain_paths,
               Catch::Matchers::RangeEquals(expected_interchain_paths));
}