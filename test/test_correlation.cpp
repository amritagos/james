#include "bondcorrel.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "fmt/core.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <optional>
#include <vector>

TEST_CASE("Test the time correlation function", "[TimeCorrelation]") {
  auto c_ij_time_series = std::vector<std::vector<int>>{
      {1, 1, 1, 1, 1}, {1, 0, 1, 1, 1}, {1, 0, 0, 1, 1},
      {1, 0, 0, 0, 1}, {1, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
  auto c_t0 =
      James::Bond::Correlation::correlation_t0_t(c_ij_time_series, 1, 1);
  INFO(
      fmt::format("Correlation at a lag time of t0, at the time origin t0={}\n",
                  c_t0.value()));
  REQUIRE_THAT(c_t0.value(), Catch::Matchers::WithinRel(1.0, 0.00001));

  // Test of the correlation function using c_ij_time_series
  std::vector<double> times{0, 10, 20, 30, 40, 50};

  // Calculate the time correlation function only for the first three elements
  // of c_ij_time_series and the times vector.
  auto c_ij_time_series_small = std::vector<std::vector<int>>{
      {1, 1, 1, 1, 1}, {1, 0, 1, 1, 1}, {1, 0, 0, 1, 1}};
  auto [tau_values, tcf_values, tcf_stderr] =
      James::Bond::Correlation::time_correlation_function(
          c_ij_time_series_small, times, 0, 1, 1, std::nullopt);
  // Required outputs
  std::vector<double> tau_expected{0, 10};
  std::vector<std::optional<double>> tcf_expected{1.0, 0.775};

  REQUIRE_THAT(tau_values, Catch::Matchers::RangeEquals(tau_expected));
  // Test that the correlation function values are as expected
  for (size_t i = 0; i < tcf_values.size(); i++) {
    REQUIRE_THAT(tcf_values[i].value(),
                 Catch::Matchers::WithinRel(tcf_expected[i].value(), 0.00001));
  }
}