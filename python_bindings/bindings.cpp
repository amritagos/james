#include "bondfinder.hpp"
#include "pairtypes.hpp"
#include "pathfinder.hpp"
#include "system.hpp"
#include "undirected_network.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

bool atom_type_found(const std::vector<int> &atom_types, int target) {
  return std::find(atom_types.begin(), atom_types.end(), target) !=
         atom_types.end();
}

std::vector<std::vector<int>>
find_interchain_ion_paths_py(const James::Atoms::System &system,
                             const std::vector<std::pair<int, int>> &pair_types,
                             std::vector<double> cutoffs,
                             const std::vector<int> &source_atom_types,
                             const std::vector<int> &destination_atom_types,
                             const std::vector<int> &intermediate_atom_types,
                             std::optional<int> max_depth,
                             bool require_intermediate, bool return_atom_ids) {

  std::vector<James::Bond::Pair> pairs{};
  pairs.reserve(pair_types.size());

  for (const auto &[type_i, type_j] : pair_types) {
    pairs.emplace_back(type_i, type_j);
  }

  auto network = Graph::UndirectedNetwork<double>(system.atoms.size());

  James::Bond::add_distance_based_bonds(network, system, pairs, cutoffs);

  auto identifier = return_atom_ids ? James::Path::WriteIdentifier::AtomID
                                    : James::Path::WriteIdentifier::Index;

  std::vector<std::vector<int>> all_paths{};

  for (size_t source_idx = 0; source_idx < system.atoms.size(); ++source_idx) {
    if (!atom_type_found(source_atom_types, system.atoms[source_idx].type)) {
      continue;
    }

    auto paths = James::Path::find_interchain_ion_pairs(
        source_idx, network, system, destination_atom_types,
        intermediate_atom_types, max_depth, require_intermediate, identifier);

    all_paths.insert(all_paths.end(), paths.begin(), paths.end());
  }

  return all_paths;
}

} // namespace

PYBIND11_MODULE(jamescpp, m) {
  m.doc() = "Python bindings for James";

  py::class_<James::Atoms::Atom>(m, "Atom")
      .def(py::init<>())
      .def_readwrite("id", &James::Atoms::Atom::id)
      .def_readwrite("type", &James::Atoms::Atom::type)
      .def_readwrite("mol_id", &James::Atoms::Atom::mol_id)
      .def_readwrite("position", &James::Atoms::Atom::position);

  py::class_<James::Atoms::System>(m, "System")
      .def(py::init<const std::vector<int> &, const std::vector<int> &,
                    const std::vector<std::vector<double>> &,
                    std::optional<std::vector<int>>,
                    std::optional<std::vector<double>>,
                    std::optional<std::vector<double>>>(),
           py::arg("ids"), py::arg("types"), py::arg("positions"),
           py::arg("mol_ids") = std::nullopt, py::arg("box") = std::nullopt,
           py::arg("box_lo") = std::nullopt)
      .def_readwrite("atoms", &James::Atoms::System::atoms)
      .def_readwrite("box", &James::Atoms::System::box)
      .def_readwrite("box_lo", &James::Atoms::System::boxLo)
      .def("n_atoms", &James::Atoms::System::n_atoms)
      .def("distance", &James::Atoms::System::distance)
      .def("collect_ids", &James::Atoms::System::collect_ids)
      .def("collect_positions", &James::Atoms::System::collect_positions);

  m.def("find_interchain_ion_paths", &find_interchain_ion_paths_py,
        py::arg("system"), py::arg("pair_types"), py::arg("cutoffs"),
        py::arg("source_atom_types"), py::arg("destination_atom_types"),
        py::arg("intermediate_atom_types"), py::arg("max_depth") = std::nullopt,
        py::arg("require_intermediate") = true,
        py::arg("return_atom_ids") = true);
}