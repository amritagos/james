from ase.data import chemical_symbols
from ase.io import read
from pathlib import Path

import solvis
from solvis.visualization import AtomicPlotter

# Input filename
script_dir = Path(__file__).resolve().parent
infilename = script_dir / "input/oct.lammpstrj"
# In the LAMMPS data file, the types of atoms are 1, 2 and 4 for O, H and Cl respectively.
fe_type = 3
cl_type = 4
h_type = 2
o_type = 1
imagename = script_dir / "temp.png"
trimmed_img_name = script_dir / "oct_no_hydrogens.png"
# Label fontsize
label_fontsize = 70

interactive_session = False
# final_camera_position = None
final_camera_position = [(24.639620235201036, 18.7699514285721, 0.3693807121895255),
 (28.208800315856934, 23.826849937438965, 21.14175033569336),
 (0.10143201079233988, -0.970477135805563, 0.21882795997142032)]

# ------------------------------------------------------
# Glean the coordinates from the LAMMPS trajectory file

# Read in the current frame
currentframe = read(
    infilename, format="lammps-dump-text"
)  # Read in the last frame of the trajectory
# Box size lengths
box_len = currentframe.get_cell_lengths_and_angles()[:3]

# Make a System object
full_system = solvis.system.System(currentframe, expand_box=True)

# Find the Fe ion (there is only one in this system)
# This will be our central atom
fe_ind = full_system.atoms.symbols.search(chemical_symbols[fe_type])
fe_pos_query_pnt = full_system.atoms.get_positions()[fe_ind[0]]

# Find solvent atoms only (exclude the Fe atom, which is the center)
solvent_atoms = full_system.atoms[
    [
        atom.index
        for atom in full_system.atoms
        if atom.number in [o_type, h_type, cl_type]
    ]
]

# Create a solvation shell, with the Fe as the center
solvation_shell = solvis.solvation_shell.create_solvation_shell_from_solvent(
    solvent_atoms, box_len, center=fe_pos_query_pnt
)
# Solvation shell without H 
solvent_atoms_without_H = full_system.atoms[
    [
        atom.index
        for atom in full_system.atoms
        if atom.number in [o_type, cl_type]
    ]
]
solvation_shell_without_H = solvis.solvation_shell.create_solvation_shell_from_solvent(
    solvent_atoms_without_H, box_len, center=fe_pos_query_pnt
)

# Create a ConvexHull object from the nearest O atom positions
convex_hull = solvation_shell.build_convex_hull_k_neighbours(
    num_neighbours=6, coordinating_type=[o_type]
)

# Calculate the sphericity from the volume and area
sph_value = solvis.util.sphericity(convex_hull.volume, convex_hull.area)
print("The calculated sphericity with 6 neighbours is", sph_value, "\n")
# ------------------------------------------------------------
# Bond and atom appearance
o_color = "midnightblue"
o_radius = 0.3
center_type = 0
central_point_colour = "black"
fe_center_radius = 0.4
cl_radius = 0.4
cl_color = "green"
h_atom_radius = 0.2
h_atom_color = "#d3d3d3"
water_bond_thickness = 0.1
special_o_color = "red"
# H bond parameters
h_segment = 0.17545664400649102
hbond_color = "orchid"
hbond_width = 40
# Decide render options for bonds
bond_opt_1 = dict(
    radius=0.15, resolution=1
)  # Fe-O bonds, single color according to bond typpe
bond_opt_2 = dict(
    bond_gradient_start=0.4, radius=water_bond_thickness, resolution=1
)  # O-H bonds, gradient, according to atomcolor type
# Render options for hydrogen bonds
hbond_opt = dict(segment_spacing=h_segment, color=hbond_color, width=hbond_width)

# Build the RendererRepresentation object
render_rep = solvis.render_helper.RendererRepresentation()
# Add the atom type and atom type specific rendering options
render_rep.add_atom_type_rendering(
    atom_type=center_type, color=central_point_colour, radius=fe_center_radius
)
render_rep.add_atom_type_rendering(atom_type=o_type, color=o_color, radius=o_radius)
render_rep.add_atom_type_rendering(
    atom_type=h_type, color=h_atom_color, radius=h_atom_radius
)
render_rep.add_atom_type_rendering(atom_type=cl_type, color=cl_color, radius=cl_radius)
# Add Fe-O bonds
render_rep.add_bond_type_rendering(bond_type=1, color=o_color, **bond_opt_1)
# Add the hydrogen bond rendering
render_rep.set_hbond_rendering(**hbond_opt)

# Add the atoms from the solvation shell
solvis.vis_initializers.fill_render_rep_atoms_from_solv_shell(
    render_rep, solvation_shell_without_H, include_center=True
)

# Add bonds from the center to the solvent atoms
bond1_render_opt = render_rep.bond_type_rendering[1].render_options
for solv_atom in solvation_shell.atoms:
    if solv_atom.number == o_type:
        iatom_tag = solv_atom.tag
        render_rep.add_bond(
            ["center", iatom_tag], bond_type=1, colorby="bondtype", **bond1_render_opt
        )

# Add the hydrogen bonds
hbond_list = solvis.bond_helper.find_all_hydrogen_bonds(
    solvation_shell,
    [o_type],
    [cl_type],
    [h_type],
    donor_H_distance=1.0,
    donor_acceptor_cutoff=3.1,
    ignore_hydrogens=True,
)
solvis.vis_initializers.fill_render_rep_hbonds_from_edges(render_rep, edges=hbond_list)

# Change the color of O atoms which are hydrogen bonded to the Cl
# easier to have hbonds with just acceptors and donors here so
hbond_list_no_h = solvis.bond_helper.find_all_hydrogen_bonds(
    solvation_shell,
    [o_type],
    [cl_type],
    [h_type],
    donor_H_distance=1.0,
    donor_acceptor_cutoff=3.1,
    ignore_hydrogens=True,
)
for bond in hbond_list_no_h:
    tag1 = bond[0]
    tag2 = bond[1]
    atom1_index = solvation_shell.tag_manager.lookup_index_by_tag(tag1)
    atom2_index = solvation_shell.tag_manager.lookup_index_by_tag(tag2)
    if solvation_shell.atoms[atom1_index].number == cl_type:
        # then update the bonded O atom
        o_name = render_rep.get_atom_name_from_tag(target_atom_tag=tag2)
        render_rep.update_atom_color(o_name, color=special_o_color)
    elif solvation_shell.atoms[atom2_index].number == cl_type:
        o_name = render_rep.get_atom_name_from_tag(target_atom_tag=tag1)
        render_rep.update_atom_color(o_name, color=special_o_color)

# Add the atom labels
# Add the labels for all solvation shell atoms, putting the indices in the currentframe as the labels 
# atom label for Fe
# Find the Fe (there is only one in this system)
fe_actor_name = render_rep.get_center_actor_name_from_type(center_type)
render_rep.set_atom_label(fe_actor_name, str(0))
# Labels for all the others (tag-1)
for atom in solvation_shell.atoms:
    if atom.number==h_type:
        continue
    index = atom.tag-1
    actor_name = render_rep.get_atom_name_from_tag(atom.tag)
    render_rep.set_atom_label(actor_name, str(index))

# Add the hull
hull_color = "#c7c7c7"
hull_options = dict(
    color=hull_color, show_edges=True, line_width=5, lighting=True, opacity=0.5
)
render_rep.add_hull(**hull_options)

# # Delete the H atoms from the solvation shell
# del solvation_shell.atoms[[atom.index for atom in solvation_shell.atoms if atom.number==h_type]]
# # Reset the atom tags
# solvation_shell.tag_manager = solvis.AtomTagManager(solvation_shell.atoms)
# ------------------------------------------------------------
# Interactive mode
if interactive_session:
    pl_inter = AtomicPlotter(interactive_mode=True, depth_peeling=True, shadows=False)

    # Populate plotter using the RendererRepresentation object
    solvis.vis_initializers.populate_plotter_from_solv_shell(
        pl_inter,
        render_rep,
        solvation_shell_without_H,
        convex_hull_list=[convex_hull],
        always_visible=True,
        show_points=False,
        font_size=label_fontsize,
    )
    # Open the interactive window
    pl_inter.interactive_window(delete_actor=True)

    # Delete the selected actors from the RendererRepresentation object
    for actor_name in pl_inter.selected_actors:
        render_rep.delete_actor(actor_name)

    final_camera_position = pl_inter.interactive_camera_position[-1]
# ------------------------------------------------------------

pl_render = AtomicPlotter(
    interactive_mode=False, window_size=[4000, 4000], depth_peeling=False, shadows=False
)
# Populate plotter using the RendererRepresentation object
solvis.vis_initializers.populate_plotter_from_solv_shell(
    pl_render,
    render_rep,
    solvation_shell_without_H,
    convex_hull_list=[convex_hull],
    always_visible=True,
    show_points=False,
    font_size=label_fontsize,
)

# Set the camera position to desired position
print("Camera position will be set to ", final_camera_position)
# Save image
pl_render.render_image(imagename, final_camera_position)

# Trim the image
trimmed_im = solvis.util.trim(imagename, 10)
trimmed_im.save(trimmed_img_name)

pl_render.close()
