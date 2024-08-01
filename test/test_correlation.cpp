#include "bondcorrel.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "fmt/core.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <vector>

TEST_CASE("Test the time correlation function", "[TimeCorrelation]") {
  auto c_ij_time_series = std::vector<std::vector<int>>{
      {1, 1, 1, 1, 1}, {1, 0, 1, 1, 1}, {1, 0, 0, 1, 1},
      {1, 0, 0, 0, 1}, {1, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
  auto c_t0 = James::Bond::Correlation::correlation_t_t0(
      c_ij_time_series, 1, 1);
  INFO(fmt::format("Correlation at a lag time of t0, at the time origin t0={}\n",
                   c_t0.value()));
  REQUIRE_THAT(c_t0.value(), Catch::Matchers::WithinRel(1.0, 0.00001));
}