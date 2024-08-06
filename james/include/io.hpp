#pragma once
#include "fmt/core.h"
#include "fstream"
#include <cstddef>
#include <fmt/ostream.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
namespace James::IO {
// Function for writing out the lag times (tau values) and the corresponding
// correlation function out to a file. zero_by_zero refers to the numeric value
// corresponding to a zero by zero TCF value.
inline void
tcf_to_file(const std::vector<double> &lag_times,
            const std::vector<std::optional<double>> &correlation_function,
            const std::string &file_path, const int zero_by_zero = 20) {
  if (lag_times.size() != correlation_function.size()) {
    throw std::invalid_argument(
        "The lag times and correlation function are not of the same size.");
  }
  double current_tcf{};
  std::fstream file;
  file.open(file_path,
            std::fstream::in | std::fstream::out | std::fstream::trunc);

  // Header string
  fmt::print(file, "# lag_time time_correlation_function\n");

  // Write the vectors out to CSV format
  for (size_t i = 0; i < lag_times.size(); i++) {
    if (!correlation_function[i].has_value()) {
      current_tcf = zero_by_zero;
    } else {
      current_tcf = correlation_function[i].has_value();
    }
    file << fmt::format("{:>25}, {:>25}\n", lag_times[i], current_tcf);
  }

  // Close the file
  file.close();
}
} // namespace James::IO