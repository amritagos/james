import pyjames as pj


def test_interchain_ion_paths_include_unphysical_motif():
    positive_residue_type = 1
    negative_residue_type = 2
    na_type = 5
    cl_type = 6

    ids = [101, 102, 201, 202, 501, 502, 503, 504, 601, 602, 301]

    types = [
        positive_residue_type,  # 0
        negative_residue_type,  # 1
        negative_residue_type,  # 2
        positive_residue_type,  # 3
        cl_type,                # 4
        na_type,                # 5
        cl_type,                # 6
        na_type,                # 7
        na_type,                # 8
        cl_type,                # 9
        negative_residue_type,  # 10
    ]

    mol_ids = [1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 3]

    positions = [
        [0.0, 0.0, 0.0],    # 0: positive residue, chain 1
        [0.0, 3.0, 0.0],    # 1: negative residue, chain 1
        [3.0, 0.0, 0.0],    # 2: negative residue, chain 2
        [9.0, 9.0, 9.0],    # 3: positive residue, chain 2, far away
        [1.0, 0.0, 0.0],    # 4: Cl
        [2.0, 0.0, 0.0],    # 5: Na
        [0.0, 1.0, 0.0],    # 6: Cl
        [0.0, 2.0, 0.0],    # 7: Na
        [-1.0, 0.0, 0.0],   # 8: Na
        [-2.0, 0.0, 0.0],   # 9: Cl
        [-3.0, 0.0, 0.0],   # 10: negative residue, chain 3
    ]

    system = pj.System(
        ids=ids,
        types=types,
        positions=positions,
        mol_ids=mol_ids,
        box=[20.0, 20.0, 20.0],
        box_lo=[0.0, 0.0, 0.0],
    )

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