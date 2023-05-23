#!/usr/bin/env python
import os
import numpy as np

from aether.diagnostics import set_diagnostics_on
from aether.indices import define_problem_type, get_ne_radius
from aether.geometry import define_node_geometry, define_1d_elements, define_rad_from_file, define_rad_from_geom
from aether.geometry import import_ply_triangles, make_data_grid
from aether.growtree import grow_tree, smooth_1d_tree
from aether.exports import export_1d_elem_geometry, export_node_geometry, export_elem_field, export_data_geometry
from aether.exports import export_1d_elem_field, export_triangle_elements, export_triangle_nodes
import pyvista as pv


########################################################################################################

def main():
    set_diagnostics_on(False)
    define_problem_type('grow_tree')  # sets up the array indices for a 1D geometry

    RUL_artery = 15  # 1d element number that supplies the RUL
    RLL_artery = 23  # 1d element number that supplies the RUL
    RML_artery = 24  # 1d element number that supplies the RUL
    LUL_artery = 7  # 1d element number that supplies the RUL
    LLL_artery = 20  # 1d element number that supplies the RUL
    angle_max = 60.0  # maximum allowed branching angle (to parent direction)
    angle_min = 20.0  # minimum allowed branching angle (to parent direction)
    branch_fraction = 0.4  # fraction of distance to COFM to branch
    length_limit = 1.0  # limit on terminal branches
    min_length = 1.2  # minimum length of elements
    rotation_limit = 180.0  # limit on angle of rotation between planes (not working)
    peel = 1.0  # proportional gap between surface and data grid (%scaling, not distance)
    # spacing = 5.48  # distance between grid points in x,y,z directions

    root_folder = '/hpc/bsha219/lung/Data'
    study = 'CTEPH'
    protocol = 'Post'
    subject = 'Alfred8'
    print("For", subject, protocol)
    input_directory = os.path.join(root_folder, study, subject, protocol, 'Vessel')
    output_directory = os.path.join(root_folder, study, subject, protocol, 'Vessel')
    ply_directory = os.path.join(root_folder, study, subject, protocol, 'Lobe', 'Lobe mesh')
    RUL_mesh = pv.read(os.path.join(ply_directory, 'RUL_surf.ply'))
    RML_mesh = pv.read(os.path.join(ply_directory, 'RML_surf.ply'))
    RLL_mesh = pv.read(os.path.join(ply_directory, 'RLL_surf.ply'))
    LUL_mesh = pv.read(os.path.join(ply_directory, 'LUL_surf.ply'))
    LLL_mesh = pv.read(os.path.join(ply_directory, 'LLL_surf.ply'))
    RUL_volume = RUL_mesh.volume
    RML_volume = RML_mesh.volume
    RLL_volume = RLL_mesh.volume
    LUL_volume = LUL_mesh.volume
    LLL_volume = LLL_mesh.volume
    # print(LLL_volume)
    RUL_spacing = (RUL_volume/6400)**(1/3)
    RML_spacing = (RML_volume/3200)**(1/3)
    RLL_spacing = (RLL_volume/8000)**(1/3)
    LUL_spacing = (LUL_volume/6400)**(1/3)
    LLL_spacing = (LLL_volume/8000)**(1/3)

    template = os.path.join(input_directory, subject + '_UpperArtery')

    define_node_geometry(template)
    define_1d_elements('/hpc/bsha219/lung/GeometricModels/3D_Digitise/Templates/art_template.ipelem')
    # define_rad_from_file(template)

    ##########################################################################
    ######################## LUL GROW ########################################
    ##########################################################################

    ply_name = 'LUL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)

    make_data_grid([0], peel, LUL_spacing, os.path.join(ply_directory, ply_name), group_name)
    export_data_geometry(os.path.join(ply_directory, ply_name), group_name, 0)

    print("Growing into LUL")

    grow_tree([0], LUL_artery, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit)
    # grow_tree([0], lingular, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit)
    # smooth_1d_tree(1, length_limit)

    ##########################################################################
    ######################## RUL GROW ########################################
    ##########################################################################

    ply_name = 'RUL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)

    make_data_grid([0], peel, RUL_spacing, os.path.join(ply_directory, ply_name), group_name)
    export_data_geometry(os.path.join(ply_directory, ply_name), group_name, 0)

    print("Growing into RUL")

    grow_tree([0], RUL_artery, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit)
    # smooth_1d_tree(1, length_limit)

    ##########################################################################
    ######################## RLL GROW ########################################
    ##########################################################################

    ply_name = 'RLL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)

    make_data_grid([0], peel, RLL_spacing, os.path.join(ply_directory, ply_name), group_name)
    export_data_geometry(os.path.join(ply_directory, ply_name), group_name, 0)

    print("Growing into RLL")
    grow_tree([0], RLL_artery, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit)
    # smooth_1d_tree(1, length_limit)

    ##########################################################################
    ######################## RML GROW ########################################
    ##########################################################################

    ply_name = 'RML_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)

    make_data_grid([0], peel, RML_spacing, os.path.join(ply_directory, ply_name), group_name)
    export_data_geometry(os.path.join(ply_directory, ply_name), group_name, 0)

    print("Growing into RML")

    grow_tree([0], RML_artery, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit)
    # smooth_1d_tree(1, length_limit)

    ##########################################################################
    ######################## LLL GROW ########################################
    ##########################################################################

    ply_name = 'LLL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)

    make_data_grid([0], peel, LLL_spacing, os.path.join(ply_directory, ply_name), group_name)
    export_data_geometry(os.path.join(ply_directory, ply_name), group_name, 0)

    print("Growing into LLL")

    grow_tree([0], LLL_artery, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit)
    # smooth_1d_tree(1, length_limit)

    ##########################################################################
    ########################  EXPORT  ########################################
    ##########################################################################

    filename = subject + '_Artery_Full'
    group_name = '1d_tree'
    export_node_geometry(os.path.join(output_directory, filename), group_name)
    export_1d_elem_geometry(os.path.join(output_directory, filename), group_name)

    # order_system = 'fit'  # fit the radii between read-in values and min_rad at order 1
    # order_options = 'all' # apply to all branches (but will only update ones with radius=0)
    # start_at = 'inlet'    # required, but previous option should make obsolete?
    # min_rad = 0.2         # radius of order 1 branches
    # h_ratio = 0.0         # doesn't matter for the 'fit' option
    # define_rad_from_geom(order_system, h_ratio, start_at, min_rad)
    # ne_radius = get_ne_radius()
    # field_name = 'radius'
    # export_1d_elem_field(ne_radius, os.path.join(output_directory, filename + '_radius'), group_name, field_name)


if __name__ == '__main__':
    main()
