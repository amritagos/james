#include "bondcorrel.hpp"
#include "io.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <filesystem>
#include <optional>
#include <vector>
namespace fs = std::filesystem;

TEST_CASE(
    "Test that you can write and read a file with the correlation function",
    "[IO]") {
  // Create a correlation function
  auto c_ij_time_series = std::vector<std::vector<int>>{
      {1, 1, 1, 1, 1}, {1, 0, 1, 1, 1}, {1, 0, 0, 1, 1},
      {1, 0, 0, 0, 1}, {1, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
  // Time series
  std::vector<double> times{0, 10, 20, 30, 40, 50};
  auto [tau_values, tcf_values, tcf_stderr] =
      James::Bond::Correlation::time_correlation_function(
          c_ij_time_series, times, 0, 1, 1, std::nullopt);

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