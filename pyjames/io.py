from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
from ase import Atoms
from ase.io import write


PathIdentifier = Literal["index", "atom_id"]


def _normalise_bonds_array(existing_bonds, n_atoms: int) -> np.ndarray:
    if existing_bonds is None:
        return np.full(n_atoms, "_", dtype=object)

    bonds = np.asarray(existing_bonds, dtype=object).copy()

    if len(bonds) != n_atoms:
        raise ValueError(
            f"Existing bonds array has length {len(bonds)}, "
            f"but atoms object has {n_atoms} atoms."
        )

    for i, entry in enumerate(bonds):
        if entry is None or str(entry).strip() == "":
            bonds[i] = "_"
        else:
            bonds[i] = str(entry)

    return bonds


def _iter_bonds(bonds: Sequence[str]):
    for atom_i, entry in enumerate(bonds):
        if entry == "_":
            continue

        for bond in str(entry).split(","):
            atom_j_str, bond_type_str = bond.split("(")
            atom_j = int(atom_j_str)
            bond_type = int(bond_type_str.rstrip(")"))
            yield atom_i, atom_j, bond_type


def _max_bond_type(bonds: Sequence[str]) -> int:
    max_type = 0

    for _, _, bond_type in _iter_bonds(bonds):
        max_type = max(max_type, bond_type)

    return max_type


def _append_bond(bonds: np.ndarray, atom_i: int, atom_j: int, bond_type: int) -> None:
    if atom_i == atom_j:
        raise ValueError("Cannot create a bond from an atom to itself.")

    # Store each undirected bond only once.
    atom_a, atom_b = sorted((atom_i, atom_j))
    new_entry = f"{atom_b}({bond_type})"

    if bonds[atom_a] == "_":
        bonds[atom_a] = new_entry
    else:
        existing_entries = str(bonds[atom_a]).split(",")
        if new_entry not in existing_entries:
            bonds[atom_a] = str(bonds[atom_a]) + "," + new_entry


def _atom_id_to_index_map(atoms: Atoms) -> dict[int, int]:
    if "id" in atoms.arrays:
        atom_ids = np.asarray(atoms.arrays["id"], dtype=int)
    else:
        atom_ids = np.arange(1, len(atoms) + 1, dtype=int)

    return {int(atom_id): index for index, atom_id in enumerate(atom_ids)}


def _path_to_indices(
    path: Sequence[int],
    atoms: Atoms,
    path_identifier: PathIdentifier,
) -> list[int]:
    if path_identifier == "index":
        indices = [int(atom_index) for atom_index in path]

    elif path_identifier == "atom_id":
        atom_id_to_index = _atom_id_to_index_map(atoms)
        indices = [atom_id_to_index[int(atom_id)] for atom_id in path]

    else:
        raise ValueError(
            f"Unknown path_identifier {path_identifier!r}. "
            "Expected 'index' or 'atom_id'."
        )

    n_atoms = len(atoms)

    for atom_index in indices:
        if atom_index < 0 or atom_index >= n_atoms:
            raise ValueError(
                f"Atom index {atom_index} is out of bounds for {n_atoms} atoms."
            )

    return indices


def write_lammps_data_with_path_bonds(
    atoms: Atoms,
    paths: Sequence[Sequence[int]],
    output_file: str | Path,
    *,
    path_identifier: PathIdentifier = "index",
    existing_bonds: Sequence[str] | None = None,
    path_bond_type: int | None = None,
    atom_style: str = "full",
    units: str = "real",
    masses: bool = True,
) -> int:
    """
    Write a LAMMPS data file with the path edges added as bonds.

    Existing bonds are kept. If existing_bonds is None, the function uses
    atoms.arrays["bonds"] if present. The path bonds are written as the next
    bond type, unless path_bond_type is given explicitly.

    Args:
        atoms:
            ASE Atoms object.
        paths:
            Paths returned by pyjames.find_interchain_ion_paths.
            Consecutive atoms in each path are written as bonds.
        output_file:
            Output LAMMPS data file.
        path_identifier:
            "index" if paths contain zero-based ASE/James atom indices.
            "atom_id" if paths contain LAMMPS atom IDs.
        existing_bonds:
            Optional ASE-style bonds array. If None, atoms.arrays["bonds"] is used
            if present.
        path_bond_type:
            Optional explicit bond type for path bonds. If None, use max existing
            bond type + 1.
        atom_style:
            LAMMPS atom style. Must be "full" for ASE to write bonds.
        units:
            LAMMPS units passed to ASE.
        masses:
            Whether ASE should write a Masses section.

    Returns:
        The bond type used for the path bonds.
    """

    if atom_style != "full":
        raise ValueError("ASE only writes bonds for atom_style='full'.")

    atoms_out = atoms.copy()
    n_atoms = len(atoms_out)

    if existing_bonds is None:
        existing_bonds = atoms.arrays.get("bonds")

    bonds = _normalise_bonds_array(existing_bonds, n_atoms)

    if path_bond_type is None:
        path_bond_type = _max_bond_type(bonds) + 1

    for path in paths:
        if len(path) < 2:
            continue

        path_indices = _path_to_indices(path, atoms_out, path_identifier)

        for atom_i, atom_j in zip(path_indices[:-1], path_indices[1:]):
            _append_bond(bonds, atom_i, atom_j, path_bond_type)

    atoms_out.arrays["bonds"] = bonds

    write(
        output_file,
        atoms_out,
        format="lammps-data",
        atom_style=atom_style,
        units=units,
        bonds=True,
        masses=masses,
    )

    return path_bond_type