#include "bondcorrel.hpp"
#include "io.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <filesystem>
#include <optional>
#include <utility>
#include <vector>
namespace fs = std::filesystem;

TEST_CASE(
    "Test that you can write and read a file with the correlation function",
    "[IO]") {
  using namespace James::Bond::Correlation;
  // Create a vector of Pairset unordered sets for each time step (5)
  PairSet set0{std::make_pair(0, 1), std::make_pair(0, 2), std::make_pair(0, 3),
               std::make_pair(1, 2), std::make_pair(1, 3)}; // At timestep 0
  PairSet set1{std::make_pair(0, 1), std::make_pair(0, 3), std::make_pair(1, 2),
               std::make_pair(1, 3)}; // At timestep 10
  PairSet set2{std::make_pair(0, 1), std::make_pair(1, 2),
               std::make_pair(1, 3)};                       // At timestep 20
  PairSet set3{std::make_pair(0, 1), std::make_pair(1, 3)}; // At timestep 30
  PairSet set4{std::make_pair(0, 1)};                       // At timestep 40
  PairSet set5{};                                           // At timestep 50
  auto pair_set_time_series =
      std::vector<PairSet>{set0, set1, set2, set3, set4, set5};
  // Time series
  std::vector<double> times{0, 10, 20, 30, 40, 50};
  auto [tau_values, tcf_values, tcf_stderr] =
      James::Bond::Correlation::time_correlation_function(
          pair_set_time_series, times, 0, 1, 1, std::nullopt);

  // Save the correlation function to a file
  auto proj_root_path = fs::current_path();
  auto tcf_file = proj_root_path / fs::path("test/tcf_out.txt");
  James::IO::tcf_to_file(tau_values, tcf_values, tcf_file, 20);

  // Read the correlation function back in
  auto [tau_values_from_file, tcf_values_from_file] =
      James::IO::tcf_from_file(tcf_file, 20);

  for (size_t i = 0; i < tau_values.size(); i++) {
    REQUIRE_THAT(tau_values_from_file[i],
                 Catch::Matchers::WithinRel(tau_values[i]));
    if (tcf_values[i].has_value()) {
      INFO(fmt::format("from file: tau = {}, tcf = {}\n",
                       tau_values_from_file[i],
                       tcf_values_from_file[i].value()));
      INFO(fmt::format("before writing: tau = {}, tcf = {}\n", tau_values[i],
                       tcf_values[i].value()));
      REQUIRE_THAT(tcf_values_from_file[i].value(),
                   Catch::Matchers::WithinRel(tcf_values[i].value()));
    } else {
      REQUIRE(tcf_values_from_file[i] == std::nullopt);
    }
  }
}