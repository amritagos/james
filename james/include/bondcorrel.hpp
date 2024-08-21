#pragma once

#include "network_base.hpp"
#include "undirected_network.hpp"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

namespace James::Bond::Correlation {
/*Creates a flattened vector containing a binary value for each pair of atoms
 * (i,j), corresponding to whether there is a connection [c(i,j)=1] or not
 * [c(i,j)=0]. The values are arranged in the order (0,1), (0,2)...(1,2),
 * (1,3)... such that i<j (equal values ommitted since the bond cannot be
 * between the atom itself).Takes the UndirectedNetwork as an
 * input. */
template <typename WeightType = double>
std::vector<int> bond_connection_info_at_tau(
    const Graph::UndirectedNetwork<WeightType> &network) {
  std::vector<int> c_ij{}; // flattened vector of c(i,j) pairs. 0 if there is no
                           // connection and 1 if there is a connection
  const size_t n_atoms = network.n_agents();
  int c_val{};

  // Loop through i j pairs
  for (size_t i = 0; i < n_atoms - 1; i++) {
    for (size_t j = i + 1; j < n_atoms; j++) {
      // For a particular pair (i,j)
      if (network.connection_exists(i, j)) {
        c_val = 1;
      } else {
        c_val = 0;
      }
      c_ij.push_back(c_val);
    }
  }

  return c_ij;
}

/*Creates a flattened vector containing a binary value for each pair of atoms
 * (i,j), corresponding to whether there is a connection [c(i,j)=1] or not
 * [c(i,j)=0]. A hydrogen bond is, however, considered to be broken if it has
 * been broken before (even if in the current step it reforms). Therefore, the
 * previous bond information is required. The values are arranged in the order
 * (0,1), (0,2)...(1,2), (1,3)... such that i<j (equal values ommitted since the
 * bond cannot be between the atom itself).Takes the UndirectedNetwork as an
 * input. */
template <typename WeightType = double>
std::vector<int> continuous_bond_connection_info_at_tau(
    const Graph::UndirectedNetwork<WeightType> &network,
    const std::vector<int> &prev_c_val) {
  std::vector<int> c_ij{}; // flattened vector of c(i,j) pairs. 0 if there is no
                           // connection and 1 if there is a connection
  const size_t n_atoms = network.n_agents();
  int c_val{static_cast<int>(prev_c_val.size())};

  // Loop through i j pairs
  size_t m = 0;
  for (size_t i = 0; i < n_atoms - 1; i++) {
    for (size_t j = i + 1; j < n_atoms; j++) {
      // For a particular pair (i,j)
      if (network.connection_exists(i, j)) {
        if (prev_c_val[m] == 0) {
          c_val = 0;
        } else {
          c_val = 1;
        }
      } else {
        c_val = 0;
      }
      c_ij.push_back(c_val);
      m++; // The position/index in the c_ij vector
    }
  }

  return c_ij;
}

/* Takes multiple `UndirectedNetwork` objects.*/
template <typename WeightType = double>
std::vector<std::vector<int>> bond_connection_info_time_series(
    const std::vector<Graph::UndirectedNetwork<WeightType>>
        &network_time_series,
    bool continuous_bond = true) {
  std::vector<std::vector<int>>
      c_ij_time_series{}; // vector of flattened vectors of c(i,j) pairs. 0 if
                          // there is no connection and 1 if there is a
                          // connection
  if (continuous_bond & (network_time_series.size() > 0)) {
    // First timestep
    auto c_ij = bond_connection_info_at_tau(network_time_series[0]);
    c_ij_time_series.push_back(c_ij);
    // Next timesteps (skip the first one)
    for (size_t i_step = 1; i_step < network_time_series.size(); i_step++) {
      c_ij = continuous_bond_connection_info_at_tau(network_time_series[i_step],
                                                    c_ij);
      c_ij_time_series.push_back(c_ij);
    }
  } else {
    for (auto &network : network_time_series) {
      c_ij_time_series.push_back(bond_connection_info_at_tau(network));
    }
  }

  return c_ij_time_series;
}

/* Calculate the correlation, at a single time origin (t0), given a lag time t
 */
template <typename T>
std::optional<double>
correlation_t0_t(const std::vector<std::vector<T>> &c_ij_time_series, size_t t0,
                 size_t t) {
  double unnorm_c_ij_val = 0.0;
  double norm_factor = 0.0;
  // Error handling
  if (t0 >= c_ij_time_series.size() || t >= c_ij_time_series.size()) {
    std::runtime_error("The time origin or lag time was greater than the "
                       "number of timesteps.\n");
  }
  if (c_ij_time_series.size() == 0) {
    std::runtime_error("Empty vector for bond information\n");
  }
  // Calculate the correlation function (average over the atoms, so over inner
  // vector size)
  for (size_t i_atom = 0; i_atom < c_ij_time_series[0].size(); i_atom++) {
    unnorm_c_ij_val +=
        c_ij_time_series[t0][i_atom] * c_ij_time_series[t][i_atom];
    norm_factor += c_ij_time_series[t0][i_atom] * c_ij_time_series[t0][i_atom];
  }
  if (norm_factor != 0) {
    return unnorm_c_ij_val / norm_factor;
  } else {
    return std::nullopt;
  }
}

/*Time correlation function, returning the tau values, the normalized
correlation function values, and the standard error in the correlation
function values.
 - start_t0: the index of the first time origin. Default value = 0
 - start_tau: first time lag (inclusive). Default value = 1.
 - delta_tau: Step size in the time lag. Default value = 1.
 - calc_upto_tau: represents the index. By default (also if set
to nullopt), this is half of the total number of steps in the c_ij_time_series
*/
template <typename T>
std::tuple<std::vector<double>, std::vector<std::optional<double>>,
           std::vector<std::optional<double>>>
time_correlation_function(const std::vector<std::vector<T>> &c_ij_time_series,
                          const std::vector<double> &time, int start_t0 = 0,
                          int start_tau = 1, int delta_tau = 1,
                          std::optional<int> calc_upto_tau = std::nullopt) {
  // Lambda to calculate the mean of a vector of TCF values (over time origins)
  // Skips std::nullopt values; returns this if the entire vector is
  // std::nullopt
  auto calc_mean = [&](const std::vector<std::optional<double>> &vec) {
    int counter = 0;   // Non 0/0 values in the vector
    double mean = 0.0; // Mean of the vector
    std::optional<double> result = std::nullopt;
    for (auto &ele : vec) {
      if (ele.has_value()) {
        counter++;
        mean += ele.value();
      }
    }
    // Now divide by the number of non 0/0 values
    if (counter > 0) {
      result = mean / counter;
    }
    // Otherwise, all values are nullopt, and that will be returned
    return result;
  };

  // Calculate the standard error
  auto calc_stderr = [&](const std::vector<std::optional<double>> &vec,
                         const std::optional<double> mean) {
    int counter = 0; // Non 0/0 values in the vector
    double std_error = 0.0;
    std::optional<double> result = std::nullopt;
    if (mean == std::nullopt) {
      return result;
    }
    for (auto &ele : vec) {
      if (ele.has_value()) {
        counter++;
        std_error +=
            (ele.value() - mean.value()) * (ele.value() - mean.value());
      }
    }
    // stderr = standard deviation/ sqrt(N)
    // corrected sample stdev = sqrt( 1/(N-1) * sum(x-x_mean)^2 )
    result = sqrt(std_error / (counter * (counter - 1)));
    return result;
  };

  // Lamba for filling up TCF values at the same lag time but different origins
  auto lag_times_over_origins =
      [&](std::vector<std::optional<double>> &tcf_values_at_lag_time,
          int current_tau) {
        for (size_t i = 0; i < tcf_values_at_lag_time.size(); i++) {
          size_t t0 = i + start_t0;
          size_t t = t0 + current_tau;
          tcf_values_at_lag_time[i] = correlation_t0_t(c_ij_time_series, t0, t);
        }
      };
  // -----------------------------
  // Handling the value of calc_upto_tau
  const int n_frames = c_ij_time_series.size();
  int max_tau = static_cast<int>(0.5 * n_frames);
  // For odd n_frames, max_tau must be inclusive, so increment by 1
  if (n_frames % 2 != 0) {
    max_tau += 1;
  }

  if (calc_upto_tau.has_value()) {
    if (calc_upto_tau.value() > max_tau) {
      std::cerr << "Warning: You set calc_upto_tau to a value greater than "
                   "max_tau. Setting to max_tau."
                << std::endl;
      calc_upto_tau = max_tau;
    }
  } else {
    calc_upto_tau = max_tau;
  }
  // -----------------------------
  // Error handling for the time
  if (time.size() < calc_upto_tau) {
    throw std::invalid_argument("The size of the time vector is incompatible "
                                "with that of the bond information array");
  }
  // -----------------------------
  std::vector<double> tau_values; // Vector of lag times
  std::vector<std::optional<double>>
      tcf_avg; // Vector of the normalized TCF, averaged over all time origins
               // for each lag time. std::nullopt when 0/0 is encountered
  std::vector<std::optional<double>>
      tcf_error; // Vector of the standard deviation errors in the TCF
  int n_origins = max_tau - start_t0; // Number of time origins
  std::vector<std::optional<double>> tcf_values_at_lag_time(
      n_origins); // Vector of the normalized TCF for a particular lag time,
                  // where each value corresponds to the lag time with respect
                  // to a time origin
  std::optional<double> mean_tcf,
      std_error_tcf; // Current value of the mean and standard error

  // Calculation of the TCF at a lag time=0
  tau_values.push_back(time[start_t0]); // update the lag time
  // Get the values of the TCF at all time origins
  lag_times_over_origins(tcf_values_at_lag_time, 0);
  // Get the mean and error, and update
  mean_tcf = calc_mean(tcf_values_at_lag_time);
  std_error_tcf = calc_stderr(tcf_values_at_lag_time, mean_tcf);
  tcf_avg.push_back(mean_tcf);
  tcf_error.push_back(std_error_tcf);

  // Calculation of the TCF at different lag times,
  // starting from start_tau
  for (int i_tau = start_tau; i_tau < calc_upto_tau.value();
       i_tau += delta_tau) {
    tau_values.push_back(time[i_tau]); // Current lag time
    // Get the TCF for this lag time over the desired time origins
    lag_times_over_origins(tcf_values_at_lag_time, i_tau);
    // Get the mean and error, and update
    mean_tcf = calc_mean(tcf_values_at_lag_time);
    std_error_tcf = calc_stderr(tcf_values_at_lag_time, mean_tcf);
    tcf_avg.push_back(mean_tcf);
    tcf_error.push_back(std_error_tcf);
  } // end of loop through lag times

  return std::make_tuple(tau_values, tcf_avg, tcf_error);
}

} // namespace James::Bond::Correlation