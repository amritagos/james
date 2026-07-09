from __future__ import annotations
from .jamescpp import Atom, System, find_interchain_ion_paths
import ase
from ase import Atoms
import numpy as np


def _get_atoms_array(atoms: Atoms, keys: list[str]):
    for key in keys:
        if key in atoms.arrays:
            return atoms.arrays[key]
    return None


def system_from_ase_atoms(atoms: Atoms) -> System:
    """Convert an ASE atoms object to a pyjames System.

    Args:
        atoms (Atoms)

    Returns:
        System
    """
    n_atoms = len(atoms)

    atom_id_array = _get_atoms_array(atoms, ["id", "atom_id", "atom-id"])

    if atom_id_array is None:
        ids = list(range(1, n_atoms + 1))
    else:
        ids = np.asarray(atom_id_array, dtype=int).tolist()

    mol_id_array = _get_atoms_array(
        atoms, ["mol-id", "mol_id", "molecule-id", "molecule_id", "mol"]
    )

    if mol_id_array is None:
        mol_ids = list(range(1, n_atoms + 1))
    else:
        mol_ids = np.asarray(mol_id_array, dtype=int).tolist()

    type_array = _get_atoms_array(atoms, ["type", "types", "atom_type", "atom_type"])

    if type_array is None:
        types = np.asarray(atoms.get_atomic_numbers(), dtype=int).tolist()
    else:
        types = np.asarray(type_array, dtype=int).tolist()

    positions = np.asarray(atoms.get_positions(), dtype=float).tolist()

    cell = atoms.cell

    if cell.rank == 3 and cell.orthorhombic:
        box = np.asarray(cell.lengths(), dtype=float).tolist()
        box_lo = [0.0, 0.0, 0.0]
    else:
        box = None
        box_lo = None

    return System(
        ids=ids,
        types=types,
        positions=positions,
        mol_ids=mol_ids,
        box=box,
        box_lo=box_lo,
    )
