#include "bondcorrel.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "fmt/core.h"
#include "undirected_network.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

using namespace James::Bond::Correlation;

Graph::UndirectedNetwork<double>
network_from_connectivity(const std::vector<size_t> &c_ij,
                          const size_t n_atoms) {
  auto network = Graph::UndirectedNetwork<double>(n_atoms);
  double weight = 1.0;
  size_t counter = 0;
  for (size_t i = 0; i < n_atoms - 1; i++) {
    for (size_t j = i + 1; j < n_atoms; j++) {
      if (c_ij[counter] == 1) {
        network.push_back_neighbour_and_weight(i, j, weight);
      }
      counter++;
    }
  }
  return network;
}

std::vector<Graph::UndirectedNetwork<double>>
create_network_time_series_test() {
  const size_t n_atoms = 4;
  std::vector<Graph::UndirectedNetwork<double>> network_time_series{};
  std::vector<std::vector<size_t>> c_ij_series{
      {1, 1, 1, 1, 1, 1}, {1, 0, 1, 1, 1, 1}, {1, 0, 0, 1, 1, 1}};

  for (const std::vector<size_t> &c_ij : c_ij_series) {
    auto network = network_from_connectivity(c_ij, n_atoms);
    network_time_series.push_back(network);
  }

  return network_time_series;
}

TEST_CASE("Test the time correlation function", "[TimeCorrelation]") {
  auto network_time_series = create_network_time_series_test();
  std::cout << "Created network time series\n";
  // Test of the correlation function using pair_set_time_series
  std::vector<double> times{0, 10, 20};

  auto [tau_values, tcf_values, tcf_stderr] =
      James::Bond::Correlation::time_correlation_function(
          network_time_series, times, 0, 1, 1, std::nullopt, true);
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