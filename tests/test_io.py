from pathlib import Path

import numpy as np
import pyjames as pj
from ase.io import read


def _parse_ase_bonds(atoms):
    edges = []

    for atom_i, entry in enumerate(atoms.arrays["bonds"]):
        if entry == "_":
            continue

        for bond in str(entry).split(","):
            atom_j_str, bond_type_str = bond.split("(")
            atom_j = int(atom_j_str)
            bond_type = int(bond_type_str.rstrip(")"))

            atom_a, atom_b = sorted((atom_i, atom_j))
            edges.append((atom_a, atom_b, bond_type))

    return sorted(edges)


def test_write_lammps_data_with_path_bonds_roundtrip(tmp_path):
    positive_residue_type = 1
    negative_residue_type = 2
    na_type = 5
    cl_type = 6

    data_file_path = Path(__file__).parent / "resources" / "test_system.data"

    atoms_ase = read(
        data_file_path,
        format="lammps-data",
        atom_style="full",
        units="real",
    )

    # Add one existing bond of type 1. ASE stores bonds as zero-based indices.
    existing_bonds = np.full(len(atoms_ase), "_", dtype=object)
    existing_bonds[0] = "1(1)"
    atoms_ase.arrays["bonds"] = existing_bonds

    system = pj.system_from_ase_atoms(atoms_ase)

    paths = pj.find_interchain_ion_paths(
        system=system,
        pair_types=[
            (positive_residue_type, cl_type),
            (negative_residue_type, na_type),
            (na_type, cl_type),
            (positive_residue_type, na_type),
            (negative_residue_type, cl_type),
        ],
        cutoffs=[1.1, 1.1, 1.1, 1.1, 1.1],
        source_atom_types=[positive_residue_type],
        destination_atom_types=[negative_residue_type],
        intermediate_atom_types=[na_type, cl_type],
        max_depth=4,
        require_intermediate=True,
        return_atom_ids=False,
    )

    assert sorted(paths) == sorted(
        [
            [2, 5, 4, 0],
            [10, 9, 8, 0],
        ]
    )

    output_file = tmp_path / "test_system_with_path_bonds.data"
    output_file = Path("test_system_with_path_bonds.data")

    path_bond_type = pj.write_lammps_data_with_path_bonds(
        atoms=atoms_ase,
        paths=paths,
        output_file=output_file,
        path_identifier="index",
        units="real",
        atom_type_count=6,
    )

    assert path_bond_type == 2
    assert output_file.exists()

    atoms_roundtrip = read(
        output_file,
        format="lammps-data",
        atom_style="full",
        units="real",
    )

    assert len(atoms_roundtrip) == 11
    assert "bonds" in atoms_roundtrip.arrays

    roundtrip_edges = _parse_ase_bonds(atoms_roundtrip)

    expected_edges = sorted(
        [
            # Existing bond, preserved as type 1.
            (0, 1, 1),
            # Path [2, 5, 4, 0], written as type 2.
            (2, 5, 2),
            (4, 5, 2),
            (0, 4, 2),
            # Path [10, 9, 8, 0], written as type 2.
            (9, 10, 2),
            (8, 9, 2),
            (0, 8, 2),
        ]
    )

    assert roundtrip_edges == expected_edges

    output_file.unlink()
    assert not output_file.exists()
