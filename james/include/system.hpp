#pragma once
#include <algorithm>
#include <optional>
#include <system_error>
#include <vector>

namespace SolvLib {

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
  std::vector<double> position{3, 0.0};     // Default initialized to 0.0
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

  // Get the number of atoms
  int n_atoms() { return atoms.size(); }
};
} // namespace SolvLib