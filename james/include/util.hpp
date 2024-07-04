#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

namespace James::Misc {
// Gets the relative distance between two positions in space
// If the box is provided, then the minimum image convention is used
inline std::vector<double> relative_distance(
    const std::vector<double> &pos1, const std::vector<double> &pos2,
    const std::optional<std::vector<double>> &box = std::nullopt) {
  auto rel_distances = std::vector<double>{};

  // Error handling
  if (pos1.size() != pos2.size()) {
    throw std::runtime_error("To calculate the relative distance, both "
                             "positions should have the same sizes");
  }
  // If the box is given, it should also be the correct size
  if (box.has_value()) {
    if (box.value().size() != pos1.size()) {
      throw std::runtime_error(
          "The box size should have the same dimensions as the position, when "
          "calculating the relative distance");
    }
  }

  const auto ndim = pos1.size();

  for (size_t i = 0; i < ndim; i++) {
    auto dr = pos1[i] - pos2[i];
    if (box.has_value()) {
      dr = dr - round(dr / box.value()[i]) * box.value()[i];
    }
    rel_distances.push_back(dr);
  }

  return rel_distances;
}
} // namespace James::Misc