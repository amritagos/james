import pyjames as pj
from pathlib import Path
from ase.io import read


def test_interchain_ion_paths_from_lammps_data_file():
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

    # ASE should preserve these LAMMPS arrays for full-style data files.
    assert "id" in atoms_ase.arrays
    assert "type" in atoms_ase.arrays
    assert "mol-id" in atoms_ase.arrays

    system = pj.system_from_ase_atoms(atoms_ase)

    assert system.n_atoms() == 11

    assert [atom.id for atom in system.atoms] == list(range(1, 12))

    assert [atom.type for atom in system.atoms] == [
        positive_residue_type,  # 1
        negative_residue_type,  # 2
        negative_residue_type,  # 3
        positive_residue_type,  # 4
        cl_type,  # 5
        na_type,  # 6
        cl_type,  # 7
        na_type,  # 8
        na_type,  # 9
        cl_type,  # 10
        negative_residue_type,  # 11
    ]

    assert [atom.mol_id for atom in system.atoms] == [1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 3]

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
        return_atom_ids=False,  # set to True for indices
    )

    # The paths contain the indices, not atom IDs
    assert sorted(paths) == sorted(
        [
            [2, 5, 4, 0],
            [10, 9, 8, 0],
        ]
    )
