#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <system_error>
#include <vector>

namespace James::Atoms {

// Holds information about an individual atom
// Assuming that the data is 3-D
struct Atom {
  Atom() = default;
  Atom(int id, int type, std::optional<int> mol_id,
       std::vector<double> &position)
      : id(id), type(type), mol_id(mol_id), position(position) {}
  int id{};   // Identifier number
  int type{}; // Type number of the atom (could be an atomic number, or type
              // number as in LAMMPS)
  std::optional<int> mol_id = std::nullopt; // Molecule identifier
  std::vector<double> position{};           // Position vector
};

class System {
public:
  std::vector<Atom> atoms{};
  std::optional<std::vector<double>> box = std::nullopt;
  std::optional<std::vector<double>> boxLo =
      std::nullopt; // Lower limits of the simulation box

  System(const std::vector<Atom> &atoms,
         std::optional<std::vector<double>> &box,
         std::optional<std::vector<double>> &boxLo)
      : atoms(atoms), box(box), boxLo(boxLo) {}

  System(const std::vector<int> &ids, const std::vector<int> &types,
         const std::vector<std::vector<double>> &positions,
         std::optional<std::vector<int>> mol_ids = std::nullopt,
         std::optional<std::vector<double>> box = std::nullopt,
         std::optional<std::vector<double>> boxLo = std::nullopt)
      : box(box), boxLo(boxLo) {
    if (mol_ids.has_value()) {
      if (ids.size() != types.size() || ids.size() != mol_ids.value().size() ||
          ids.size() != positions.size()) {
        throw std::invalid_argument("Input vectors must have the same size");
      }
      for (size_t i = 0; i < ids.size(); ++i) {
        atoms.emplace_back(ids[i], types[i], mol_ids.value()[i],
                           const_cast<std::vector<double> &>(positions[i]));
      }
    } else {
      if (ids.size() != types.size() || ids.size() != positions.size()) {
        throw std::invalid_argument("Input vectors must have the same size");
      }
      for (size_t i = 0; i < ids.size(); ++i) {
        atoms.emplace_back(ids[i], types[i], std::nullopt,
                           const_cast<std::vector<double> &>(positions[i]));
      }
    }
  }

  System() = default;

  // Delete the (n+1)^th Atom in the System object
  void del(int n) { atoms.erase(atoms.begin() + n); }

  // Delete a range of Atom objects, in the range [first, last)
  // For instance, if first=0 and last=3, the first three Atom objects will be
  // deleted
  void del(int first, int last) {
    atoms.erase(atoms.begin() + first, atoms.begin() + last);
  }

  // Add an atom to the System object
  void push_back(const Atom &atom) { atoms.push_back(atom); }

  // Collect all the IDs in the System object
  std::vector<int> collect_ids() {
    auto ids = std::vector<int>{};

    std::transform(atoms.begin(), atoms.end(), std::back_inserter(ids),
                   [](Atom const &a) -> double { return a.id; });

    return ids;
  }

  // Get the index corresponding to the atom ID
  std::optional<size_t> index_from_id(int target_id) const {

    auto is_target_id = [&target_id](Atom atom) {
      if (atom.id == target_id) {
        return true;
      } else {
        return false;
      }
    };

    // Search for the target ID
    if (auto it = std::find_if(atoms.begin(), atoms.end(), is_target_id);
        it != atoms.end())
      return std::distance(atoms.begin(), it);
    else
      return std::nullopt;
  }

  // Find all the indices (not IDs actually) in atoms with the desired molecule
  // ID
  std::vector<size_t> find_atoms_in_molecule(int target_mol_id) const {
    std::vector<size_t> indices{};
    auto it = atoms.begin();

    while (it != atoms.end()) {
      it = std::find_if(it, atoms.end(), [target_mol_id](const Atom &atom) {
        if (atom.mol_id.has_value()) {
          return atom.mol_id == target_mol_id;
        } else {
          return false;
        }
      });

      if (it != atoms.end()) {
        indices.push_back(std::distance(atoms.begin(), it));
        ++it; // Move iterator to the next element
      }
    }

    return indices;
  }

  // Distance (with or without periodic boundary conditions) between two Atom
  // objects in the System
  double distance(size_t index1, size_t index2) const {
    double r = 0.0;
    const size_t ndim = 3; // Number of dimensions (3)
    // For three dimensions
    for (size_t i = 0; i < ndim; i++) {
      auto dr = atoms[index1].position[i] - atoms[index2].position[i];
      if (box.has_value()) {
        dr = dr - round(dr / box.value()[i]) * box.value()[i];
      }
      r += dr * dr; // Update the distance
    }

    return sqrt(r);
  }

  // Get the number of atoms
  int n_atoms() { return atoms.size(); }
};
} // namespace James::Atoms