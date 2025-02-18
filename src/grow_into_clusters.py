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

    RUL_artery = 15  # 1d element number that supplies the RUL
    RLL_artery = 23  # 1d element number that supplies the RLL
    RML_artery = 24  # 1d element number that supplies the RML
    LUL_artery = [11, 19]  # 1d element number that supplies the LUL
    LLL_artery = 20  # 1d element number that supplies the LLL
    angle_max = 60.0  # maximum allowed branching angle (to parent direction)
    angle_min = 20.0  # minimum allowed branching angle (to parent direction)
    branch_fraction = 0.4  # fraction of distance to COFM to branch
    length_limit = 1.0  # limit on terminal branches
    min_length = 1.2  # minimum length of elements
    rotation_limit = 180.0  # limit on angle of rotation between planes (not working)
    peel = 1.0  # proportional gap between surface and data grid (%scaling, not distance)
    # spacing = 5.48  # distance between grid points in x,y,z directions

    # Mapping clusters to elems in lobes
    branch_mapping = {
        1: {'Lobe': 'RUL', 'branch_name': 'Posterior_lower', 'element': 52},
        2: {'Lobe': 'RUL', 'branch_name': 'Apical_upper', 'element': 50},
        3: {'Lobe': 'RUL', 'branch_name': 'Apical_lower', 'element': 49},
        4: {'Lobe': 'RUL', 'branch_name': 'Anterior_lower', 'element': 35},
        5: {'Lobe': 'RUL', 'branch_name': 'Anterior_upper', 'element': 36},
        6: {'Lobe': 'RUL', 'branch_name': 'Posterior_upper', 'element': 51},
        7: {'Lobe': 'RML', 'branch_name': 'Medial_posterior', 'element': 60},
        8: {'Lobe': 'RML', 'branch_name': 'Lateral_anterior', 'element': 58},
        9: {'Lobe': 'RML', 'branch_name': 'Medial_anterior', 'element': 59},
        10: {'Lobe': 'RML', 'branch_name': 'Lateral_posterior', 'element': 57},
        11: {'Lobe': 'RLL', 'branch_name': 'Posterior', 'element': 67},
        12: {'Lobe': 'RLL', 'branch_name': 'Superior_upper', 'element': 54},
        13: {'Lobe': 'RLL', 'branch_name': 'Anterior', 'element': 65},
        14: {'Lobe': 'RLL', 'branch_name': 'Superior_lower', 'element': 53},
        15: {'Lobe': 'RLL', 'branch_name': 'Medial', 'element': 68},
        16: {'Lobe': 'RLL', 'branch_name': 'Lateral', 'element': 66},
        17: {'Lobe': 'LLL', 'branch_name': 'Posterior', 'element': 61},
        18: {'Lobe': 'LLL', 'branch_name': 'Anterior', 'element': 46},
        19: {'Lobe': 'LLL', 'branch_name': 'Superior_upper', 'element': 63},
        20: {'Lobe': 'LLL', 'branch_name': 'Lateral', 'element': 45},
        21: {'Lobe': 'LLL', 'branch_name': 'Superior_lower', 'element': 64},
        22: {'Lobe': 'LLL', 'branch_name': 'Medial', 'element': 62},
        23: {'Lobe': 'LUL', 'branch_name': 'Anterior_upper', 'element': 27},
        24: {'Lobe': 'LUL', 'branch_name': 'Apical_lower', 'element': 41},
        25: {'Lobe': 'LUL', 'branch_name': 'Posterior_lower', 'element': 44},
        26: {'Lobe': 'LUL', 'branch_name': 'Anterior_upper', 'element': 28},
        27: {'Lobe': 'LUL', 'branch_name': 'Posterior_upper', 'element': 43},
        28: {'Lobe': 'LUL', 'branch_name': 'Lingular_upper', 'element': 29},
        29: {'Lobe': 'LUL', 'branch_name': 'Lingular_lower', 'element': 30},
        30: {'Lobe': 'LUL', 'branch_name': 'Apical_upper', 'element': 42},
    }

    root_folder = '/hpc/bsha219/lung/Data'
    study = 'CTEPH'
    protocol = 'Pre'
    subject = 'Alfred1'
    print("For", subject, protocol)
    input_directory = os.path.join(root_folder, study, subject, protocol, 'Vessel')
    output_directory = os.path.join(root_folder, study, subject, protocol, 'Vessel')
    lobe_directory = os.path.join(root_folder, study, subject, protocol, 'Lobe')
    intensity_path = os.path.join(root_folder, study, subject, protocol, 'Intensity_mapping')
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

    template = os.path.join(input_directory, subject + '_UpperArtery')

    define_node_geometry(template)
    define_1d_elements('/hpc/bsha219/lung/GeometricModels/3D_Digitise/Templates/art_template.ipelem')
    # define_rad_from_file(template)

    ##########################################################################
    ######################## LUL GROW ########################################
    ##########################################################################

    # Loop over each cluster in the mapping and call the two functions.
    for cluster, info in branch_mapping.items():
        lobe = info['Lobe']
        branch_name = info['branch_name']
        element = info['element']
        print("Growing {} {}".format(lobe, branch_name))
        ply_name = str(lobe) + '_surf'
        group_name = ply_name
        import_ply_triangles(os.path.join(ply_directory, ply_name))
        export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
        export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)

        # Construct the filename for define_data_geometry:
        # This will yield something like 'Cl_1.ipdata', 'Cl_2.ipdata', etc.
        geometry_filename = f"Cl_{cluster}.ipdata"
        geometry_filepath = os.path.join(intensity_path, geometry_filename)

        # Construct the mapping filename for grow_tree:
        # This will yield something like 'RUL_posterior_lower_mapping'
        # (Note: converting branch_name to lowercase to match your example.)
        mapping_filename = f"{lobe}_{branch_name.lower()}_mapping"
        mapping_filepath = os.path.join(lobe_directory, mapping_filename)

        # Call the first function with the geometry file path.
        define_data_geometry(geometry_filepath)

        # Call the second function with the appropriate parameters.
        # Here, we use the element number in place of the second argument.
        # print(element)
        # print(cluster)
        # print(type(element))
        grow_tree([0], element, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
                  True, mapping_filepath, 'closest')

    # lobe = 'LUL'
    # ply_name = 'LUL_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(ply_directory, ply_name))
    # export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(lobe_directory, ply_name), group_name)
    # print("Growing into LUL")
    # # Anterior Lower
    # define_data_geometry(os.path.join(intensity_path, 'Cl_28.ipdata'))
    # grow_tree([0], 28, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LUL_anterior_lower_mapping'), 'closest')
    # # Anterior upper
    # define_data_geometry(os.path.join(intensity_path, 'Cl_30.ipdata'))
    # grow_tree([0], 27, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LUL_anterior_upper_mapping'), 'closest')
    # # Posterior upper
    # define_data_geometry(os.path.join(intensity_path, 'Cl_25.ipdata'))
    # grow_tree([0], 43, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LUL_posterior_upper_mapping'), 'closest')
    # # Posterior lower
    # define_data_geometry(os.path.join(intensity_path, 'Cl_26.ipdata'))
    # grow_tree([0], 44, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LUL_posterior_lower_mapping'), 'closest')
    # # Apical upper
    # define_data_geometry(os.path.join(intensity_path, 'Cl_23.ipdata'))
    # grow_tree([0], 41, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LUL_apical_upper_mapping'), 'closest')
    # # Apical lower
    # define_data_geometry(os.path.join(intensity_path, 'Cl_29.ipdata'))
    # grow_tree([0], 42, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LUL_apical_lower_mapping'), 'closest')
    # # Lingular lateral
    # define_data_geometry(os.path.join(intensity_path, 'Cl_24.ipdata'))
    # grow_tree([0], 29, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LUL_lingular_medial_mapping'), 'closest')
    # # Lingular medial
    # define_data_geometry(os.path.join(intensity_path, 'Cl_27.ipdata'))
    # grow_tree([0], 30, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LUL_lingular_lateral_mapping'), 'closest')
    #
    # ##########################################################################
    # ######################## RUL GROW ########################################
    # ##########################################################################
    #
    # lobe = 'RUL'
    # ply_name = 'RUL_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(ply_directory, ply_name))
    # export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)
    # print("Growing into RUL")
    # # Apical upper
    # define_data_geometry(os.path.join(intensity_path, 'Cl_5.ipdata'))
    # grow_tree([0], 49, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RUL_apical_upper_mapping'), 'closest')
    # # Apical lower
    # define_data_geometry(os.path.join(intensity_path, 'Cl_3.ipdata'))
    # grow_tree([0], 50, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RUL_apical_lower_mapping'), 'closest')
    # # anterior upper
    # define_data_geometry(os.path.join(intensity_path, 'Cl_6.ipdata'))
    # grow_tree([0], 36, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RUL_anterior_lower_mapping'), 'closest')
    # # Anterior lower
    # define_data_geometry(os.path.join(intensity_path, 'Cl_2.ipdata'))
    # grow_tree([0], 35, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RUL_anterior_upper_mapping'), 'closest')
    # # Posterior upper
    # define_data_geometry(os.path.join(intensity_path, 'Cl_1.ipdata'))
    # grow_tree([0], 51, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RUL_posterior_upper_mapping'), 'closest')
    # # Posterior lower
    # define_data_geometry(os.path.join(intensity_path, 'Cl_4.ipdata'))
    # grow_tree([0], 52, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RUL_posterior_lower_mapping'), 'closest')
    #
    # ##########################################################################
    # ######################## RLL GROW ########################################
    # ##########################################################################
    #
    # lobe = 'RLL'
    # ply_name = 'RLL_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(ply_directory, ply_name))
    # export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)
    # print("Growing into RLL")
    # # Medial
    # define_data_geometry(os.path.join(intensity_path, 'Cl_15.ipdata'))
    # grow_tree([0], 65, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RLL_lateral_mapping'), 'closest')
    # # Lateral
    # define_data_geometry(os.path.join(intensity_path, 'Cl_14.ipdata'))
    # grow_tree([0], 67, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RLL_anterior_mapping'), 'closest')
    # # Anterior
    # define_data_geometry(os.path.join(intensity_path, 'Cl_16.ipdata'))
    # grow_tree([0], 66, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RLL_posterior_mapping'), 'closest')
    # # Posterior
    # define_data_geometry(os.path.join(intensity_path, 'Cl_12.ipdata'))
    # grow_tree([0], 68, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RLL_medial_mapping'), 'closest')
    # # Superior upper
    # define_data_geometry(os.path.join(intensity_path, 'Cl_13.ipdata'))
    # grow_tree([0], 54, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RLL_superior_lower_mapping'), 'closest')
    # # Superior lower
    # define_data_geometry(os.path.join(intensity_path, 'Cl_11.ipdata'))
    # grow_tree([0], 53, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RLL_superior_upper_mapping'), 'closest')
    #
    # ##########################################################################
    # ######################## RML GROW ########################################
    # ##########################################################################
    #
    # lobe = 'RML'
    # ply_name = 'RML_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(ply_directory, ply_name))
    # export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)
    # print("Growing into RML")
    # # Medial anterior
    # define_data_geometry(os.path.join(intensity_path, 'Cl_8.ipdata'))
    # grow_tree([0], 60, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RML_medial_anterior_mapping'), 'closest')
    # # Medial posterior
    # define_data_geometry(os.path.join(intensity_path, 'Cl_7.ipdata'))
    # grow_tree([0], 59, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RML_medial_posterior_mapping'), 'closest')
    # # Lateral Anterior
    # define_data_geometry(os.path.join(intensity_path, 'Cl_9.ipdata'))
    # grow_tree([0], 58, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RML_lateral_posterior_mapping'), 'closest')
    # # Lateral Posterior
    # define_data_geometry(os.path.join(intensity_path, 'Cl_10.ipdata'))
    # grow_tree([0], 57, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'RML_lateral_anterior_mapping'), 'closest')
    #
    # ##########################################################################
    # ######################## LLL GROW ########################################
    # ##########################################################################
    #
    # lobe = 'LLL'
    # ply_name = 'LLL_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(ply_directory, ply_name))
    # export_triangle_nodes(os.path.join(ply_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(ply_directory, ply_name), group_name)
    # print("Growing into LLL")
    # # Medial
    # define_data_geometry(os.path.join(intensity_path, 'Cl_17.ipdata'))
    # grow_tree([0], 62, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LLL_lateral_mapping'), 'closest')
    # # Lateral
    # define_data_geometry(os.path.join(intensity_path, 'Cl_21.ipdata'))
    # grow_tree([0], 45, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LLL_superior_lower_mapping'), 'closest')
    # # Anterior
    # define_data_geometry(os.path.join(intensity_path, 'Cl_19.ipdata'))
    # grow_tree([0], 46, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LLL_superior_upper_mapping'), 'closest')
    # # Posterior
    # define_data_geometry(os.path.join(intensity_path, 'Cl_18.ipdata'))
    # grow_tree([0], 61, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LLL_anterior_mapping'), 'closest')
    # # Superior upper
    # define_data_geometry(os.path.join(intensity_path, 'Cl_22.ipdata'))
    # grow_tree([0], 64, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LLL_posterior_mapping'), 'closest')
    # # Superior lower
    # define_data_geometry(os.path.join(intensity_path, 'Cl_20.ipdata'))
    # grow_tree([0], 63, 0, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit,
    #           True, os.path.join(lobe_directory, 'LLL_medial_mapping'), 'closest')

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
