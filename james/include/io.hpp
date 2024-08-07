#pragma once
#include "fmt/core.h"
#include "fstream"
#include <cstddef>
#include <filesystem>
#include <fmt/ostream.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace James::IO {
// Read the entire file into std::string
// Borrowed from seldon
inline std::string get_file_contents(const std::string &filename) {
  auto path = std::filesystem::path(filename);
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error(fmt::format(
        "Canot read from {}. File does not exist!", fmt::streamed(path)));
  }

  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in) {
    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return (contents);
  }
  throw(std::runtime_error("Cannot read file."));
}
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
  std::fstream file;
  file.open(file_path,
            std::fstream::in | std::fstream::out | std::fstream::trunc);

  // Header string
  fmt::print(file, "# lag_time time_correlation_function\n");

  // Write the vectors out to CSV format
  for (size_t i = 0; i < lag_times.size(); i++) {
    if (correlation_function[i].has_value()) {
      file << fmt::format("{:>25}, {:>25}\n", lag_times[i],
                          correlation_function[i].value());
    } else {
      file << fmt::format("{:>25}, {:>25}\n", lag_times[i], zero_by_zero);
    }
  }

  // Close the file
  file.close();
}
inline std::pair<std::vector<double>, std::vector<std::optional<double>>>
tcf_from_file(const std::string &file_path, const int zero_by_zero = 20) {
  std::vector<double> lag_times{};
  std::vector<std::optional<double>> correlation_function{};

  // Get the entire file contents as a string
  std::string file_contents = get_file_contents(file_path);

  // bool finished = false;
  size_t start_of_line = 0;
  bool finished = false;
  while (!finished) {
    // Find the end of the current line
    auto end_of_line = file_contents.find('\n', start_of_line);
    if (end_of_line == std::string::npos) {
      finished = true;
    }

    // Get the current line as a substring
    auto line =
        file_contents.substr(start_of_line, end_of_line - start_of_line);
    start_of_line = end_of_line + 1;

    if (line.empty()) {
      break;
    }

    if (line[0] == '#') {
      continue;
    }

    // Parse the columns

    //@TODO: refactor with util/parse_comma_separated_list
    size_t start_of_column = 0;
    bool finished_row = false;
    size_t idx_column = 0;
    while (!finished_row) {
      auto end_of_column = line.find(',', start_of_column);

      if (end_of_column == std::string::npos) {
        finished_row = true;
      }

      auto column_substring =
          line.substr(start_of_column, end_of_column - start_of_column);
      start_of_column = end_of_column + 1;

      // First column is the lag time
      if (idx_column == 0) {
        lag_times.push_back(std::stod(column_substring));
      }
      // The second column contains the correlation function
      else if (idx_column == 1) {
        const auto tcf_value = std::stod(column_substring);
        // If tcf_value is the zero_by_zero value then it is nullopt
        if (tcf_value == zero_by_zero) {
          correlation_function.push_back(std::nullopt);
        } else {
          correlation_function.push_back(tcf_value);
        }
      }
      idx_column++;
    }
  }

  // Return the lag times and correlation function
  return std::make_pair(lag_times, correlation_function);
}
} // namespace James::IO