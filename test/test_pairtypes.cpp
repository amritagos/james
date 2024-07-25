#include "pairtypes.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <vector>

TEST_CASE("Test the Pair class", "[Pair]") {
  auto pair31 = James::Bond::Pair(3, 1);
  auto pair13 = James::Bond::Pair(1, 3);

  // Commutative pairs should be the same
  REQUIRE(pair13 == pair31);
}

TEST_CASE("Test that you can make an unordered_map with Pairs as keys",
          "[PairCutoffMap]") {
  auto pair13 = James::Bond::Pair(1, 3);
  auto pair11 = James::Bond::Pair(1, 1);
  auto pairs = std::vector<James::Bond::Pair>{pair13, pair11};
  auto cutoffs = std::vector<double>{3.2, 4.0};
  auto pair_cutoff_map = James::Bond::create_pairtype_cutoffs(pairs, cutoffs);

  REQUIRE(pair_cutoff_map[pair13] == cutoffs[0]);
  REQUIRE(pair_cutoff_map[pair11] == cutoffs[1]);

  // Update 1-3
  pair_cutoff_map[James::Bond::Pair(3, 1)] = 3.5;
  REQUIRE(pair_cutoff_map[pair13] == 3.5);
}