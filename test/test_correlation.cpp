#include "bondcorrel.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "fmt/core.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

using namespace James::Bond::Correlation;
std::vector<PairSet> create_pair_set_time_series_test() {
  const size_t n_atoms = 4;
  std::vector<PairSet> pair_set_time_series{};
  PairSet pair_set{};

  // All possible pairs (4 atoms here)
  for (size_t i = 0; i < n_atoms - 1; i++) {
    for (size_t j = i + 1; j < n_atoms; j++) {
      pair_set.emplace(i, j);
    }
  }

  // First timestep {1, 1, 1, 1, 1, 1}
  pair_set_time_series.push_back(pair_set);

  // Second timestep {1, 0, 1, 1, 1, 1}
  pair_set.erase(std::make_pair(0, 2));
  pair_set_time_series.push_back(pair_set);

  // Third timestep {1, 0, 0, 1, 1, 1}
  pair_set.erase(std::make_pair(0, 3));
  pair_set_time_series.push_back(pair_set);

  return pair_set_time_series;
}

TEST_CASE("Test the time correlation function", "[TimeCorrelation]") {
  auto pair_set_time_series = create_pair_set_time_series_test();
  REQUIRE(pair_set_time_series.size() == 3);
  auto c_t0 =
      James::Bond::Correlation::correlation_t0_t(pair_set_time_series, 1, 1);
  INFO(
      fmt::format("Correlation at a lag time of t0, at the time origin t0={}\n",
                  c_t0.value()));
  REQUIRE_THAT(c_t0.value(), Catch::Matchers::WithinRel(1.0, 0.00001));

  // Test of the correlation function using pair_set_time_series
  std::vector<double> times{0, 10, 20};

  auto [tau_values, tcf_values, tcf_stderr] =
      James::Bond::Correlation::time_correlation_function(
          pair_set_time_series, times, 0, 1, 1, std::nullopt);
  // Required outputs
  std::vector<double> tau_expected{0, 10};
  std::vector<std::optional<double>> tcf_expected{1.0, 0.81666667};

  REQUIRE_THAT(tau_values, Catch::Matchers::RangeEquals(tau_expected));
  // Test that the correlation function values are as expected
  for (size_t i = 0; i < tcf_expected.size(); i++) {
    REQUIRE_THAT(tcf_values[i].value(),
                 Catch::Matchers::WithinRel(tcf_expected[i].value(), 0.00001));
  }
}