#pragma once

#include "fmt/core.h"
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

// Function to calculate the angle ABC (in radians) given the coordinates for
// points A, B and C.
inline double
angleABC(const std::vector<double> &A, const std::vector<double> &B,
         const std::vector<double> &C,
         const std::optional<std::vector<double>> &box = std::nullopt) {
  // Ensure the input vectors have exactly the same number of elements
  if (A.size() != B.size() || B.size() != C.size()) {
    throw std::invalid_argument(
        "All points must have exactly the same number of dimensions.");
  }

  // Vector AB
  auto ab_vector = relative_distance(A, B, box);

  // Vector CB
  auto cb_vector = relative_distance(C, B, box);

  // Dot product of AB and CB
  double dot_product = 0.0;
  // Magnitudes of AB and CB
  double ab_norm = 0.0;
  double cb_norm = 0.0;
  for (size_t i = 0; i < ab_vector.size(); i++) {
    dot_product += ab_vector[i] * cb_vector[i];
    ab_norm += ab_vector[i] * ab_vector[i];
    cb_norm += cb_vector[i] * cb_vector[i];
  }
  ab_norm = std::sqrt(ab_norm);
  cb_norm = std::sqrt(cb_norm);

  // Check for zero magnitude to avoid division by zero
  if (ab_norm == 0.0 || cb_norm == 0.0) {
    throw std::runtime_error("One of the vectors has zero length.");
  }

  // Cosine of the angle
  double cos_theta = dot_product / (ab_norm * cb_norm);

  // To avoid any numerical errors leading to values slightly outside [-1, 1]
  if (cos_theta < -1.0)
    cos_theta = -1.0;
  if (cos_theta > 1.0)
    cos_theta = 1.0;

  // Angle in radians
  return std::acos(cos_theta);
}

// Function to calculate the angle ABC (in degrees) given the coordinates for
// points A, B and C.
inline double
angleABCdeg(const std::vector<double> &A, const std::vector<double> &B,
            const std::vector<double> &C,
            const std::optional<std::vector<double>> &box = std::nullopt) {
  // Convert to degrees
  double angleDegrees = angleABC(A, B, C, box) * (180.0 / M_PI);
  return angleDegrees;
}

// Downsample a vector such that skip_every steps are skipped. Can only accept
// values greater than 0.
template <typename T>
std::vector<T> downsample(const std::vector<T> &vec, int skip_every) {
  // Create a new vector to hold the downsampled vector
  auto results = std::vector<T>{};
  // Error handling of skip_every
  if (skip_every <= 0) {
    throw std::invalid_argument(
        fmt::format("You entered an invalid value of {} for skip_every, which "
                    "must be 1 or greater.",
                    skip_every));
  }

  for (size_t i = 0; i < vec.size(); i++) {
    results.push_back(vec[i]);
  }

  return results;
}

} // namespace James::Misc