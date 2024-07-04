#include "system.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
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
  system.push_back(James::Atoms::Atom(2, 4, 2, positions[1]));

  // Test that the distance between the O and Cl is 3.02442
  double dist_required = 3.02442;

  REQUIRE_THAT(system.distance(0, 1),
               Catch::Matchers::WithinRel(dist_required, 1e-5));
}