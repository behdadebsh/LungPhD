#!/usr/bin/env python
import os
import numpy as np

from aether.diagnostics import set_diagnostics_on
from aether.indices import define_problem_type, get_ne_radius
from aether.geometry import define_node_geometry, define_1d_elements, define_rad_from_file, define_rad_from_geom, define_data_geometry
from aether.geometry import import_ply_triangles, make_data_grid
from aether.growtree import grow_tree, smooth_1d_tree
from aether.exports import export_1d_elem_geometry, export_node_geometry, export_elem_field, export_data_geometry
from aether.exports import export_1d_elem_field, export_triangle_elements, export_triangle_nodes
import pyvista as pv


########################################################################################################

def main():
    set_diagnostics_on(False)
    define_problem_type('grow_tree')  # sets up the array indices for a 1D geometry

    peel = 1.0  # proportional gap between surface and data grid (%scaling, not distance)

    root_folder = '/hpc/bsha219/lung/Data'
    study = 'CTEPH'
    protocol = 'Pre'
    subject = 'Alfred5'
    print("For", subject, protocol)
    lobe_directory = os.path.join(root_folder, study, subject, protocol, 'Lobe')
    ply_directory = os.path.join(root_folder, study, subject, protocol, 'Lobe', 'Lobe_mesh')
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

    RUL_spacing = (RUL_volume/6400)**(1/3)
    RML_spacing = (RML_volume/3200)**(1/3)
    RLL_spacing = (RLL_volume/8000)**(1/3)
    LUL_spacing = (LUL_volume/6400)**(1/3)
    LLL_spacing = (LLL_volume/8000)**(1/3)

    ##########################################################################
    ######################## LUL  ########################################
    ##########################################################################

    lobe = 'LUL'
    ply_name = 'LUL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(lobe_directory, ply_name), group_name)
    print(f'For {lobe}')
    make_data_grid([0], 0, peel, LUL_spacing)
    export_data_geometry(os.path.join(lobe_directory, lobe + '_datapoints'), group_name, 0)

    ##########################################################################
    ######################## RUL ########################################
    ##########################################################################

    lobe = 'RUL'
    ply_name = 'RUL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)
    print(f'For {lobe}')
    make_data_grid([0], 0, peel, RUL_spacing)
    export_data_geometry(os.path.join(lobe_directory, lobe + '_datapoints'), group_name, 0)

    ##########################################################################
    ######################## RLL ########################################
    ##########################################################################

    lobe = 'RLL'
    ply_name = 'RLL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)
    print(f'For {lobe}')
    make_data_grid([0], 0, peel, RLL_spacing)
    export_data_geometry(os.path.join(lobe_directory, lobe + '_datapoints'), group_name, 0)

    ##########################################################################
    ######################## RML ########################################
    ##########################################################################

    lobe = 'RML'
    ply_name = 'RML_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)
    print(f'For {lobe}')
    make_data_grid([0], 0, peel, RML_spacing)
    export_data_geometry(os.path.join(lobe_directory, lobe + '_datapoints'), group_name, 0)

    ##########################################################################
    ######################## LLL ########################################
    ##########################################################################

    lobe = 'LLL'
    ply_name = 'LLL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(ply_directory, ply_name))
    export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)
    print(f'For {lobe}')
    make_data_grid([0], 0, peel, LLL_spacing)
    export_data_geometry(os.path.join(lobe_directory, lobe + '_datapoints'), group_name, 0)


if __name__ == '__main__':
    main()
