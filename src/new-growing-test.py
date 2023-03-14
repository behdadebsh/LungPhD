from aether.diagnostics import set_diagnostics_on
from aether.indices import define_problem_type
from aether.geometry import define_node_geometry_2d, define_elem_geometry_2d, make_data_grid, evaluate_ordering, \
    define_data_geometry, define_node_geometry, define_1d_elements, group_elem_parent_term
from aether.exports import export_node_geometry_2d, export_elem_geometry_2d, export_data_geometry, export_node_geometry, \
    export_1d_elem_geometry
from aether.growtree import grow_tree
import os


set_diagnostics_on(False)
define_problem_type('grow_tree')

study = 'CTEPH'
subject = 'CTEPH10'
volume = 'FRC'
subject_path = os.path.join('/hpc/bsha219/lung/Data/', study, subject, volume)
apply_mapping = False   # default should be false for old growing - Make True for new grow to clusters
# subject_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH3/FRC'

define_node_geometry(os.path.join(subject_path, 'Vessel', subject + '_UpperArtery.ipnode'))
# define_node_geometry('/hpc/bsha219/lung/Data/CTEPH/NewCase/FRC/Vessel/NewCase_UpperArtery.ipnode')
define_1d_elements('/hpc/bsha219/lung/GeometricModels/3D_Digitise/Templates/art_template.ipelem')

print('========= 1D mesh read ========')

# Growing into RLL
print("=========== Growing into RLL ==========")

define_node_geometry_2d(os.path.join(subject_path, 'Lung', 'SurfaceFEMesh', 'Right_fitted.ipnode'))
# define_node_geometry_2d('/hpc/bsha219/lung/Data/CTEPH/NewCase/FRC/Lung/SurfaceFEMesh/Right_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/RLL_surface.ipelem', 'unit')
print("=========== RLL mesh read ==========")
if apply_mapping:
    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping','RLL_CL1.ipdata'))
    print("============ RLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(53)
    print("======== Grouped by parent element #53 =========")
    grow_tree(53, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RLL_CL2_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RLL_CL6.ipdata'))
    print("============ RLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(54)
    print("======== Grouped by parent element #54 =========")
    grow_tree(54, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RLL_CL4_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RLL_CL4.ipdata'))
    print("============ RLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(65)
    print("======== Grouped by parent element #65 =========")
    grow_tree(65, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RLL_CL1_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RLL_CL3.ipdata'))
    print("============ RLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(66)
    print("======== Grouped by parent element #66 =========")
    grow_tree(66, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RLL_CL4_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RLL_CL5.ipdata'))
    print("============ RLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(67)
    print("======== Grouped by parent element #67 =========")
    grow_tree(67, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RLL_CL5_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RLL_CL2.ipdata'))
    print("============ RLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(68)
    print("======== Grouped by parent element #68 =========")
    grow_tree(68, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RLL_CL3_mapping')
    print("==== RLL Grown ======")
else:
    make_data_grid(0, 4.5, False, 'data_grid', 'data_grid')
    print("============ RLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(23)
    print("======== Grouped by parent element #23 =========")
    grow_tree(23, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, True, 'RLL_mapping')
    print("==== RLL Grown ======")


# Growing into RML
print("=========== Growing into RML ==========")

define_node_geometry_2d(os.path.join(subject_path, 'Lung', 'SurfaceFEMesh', 'Right_fitted.ipnode'))
# define_node_geometry_2d('/hpc/bsha219/lung/Data/CTEPH/NewCase/FRC/Lung/SurfaceFEMesh/Right_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/RML_surface.ipelem', 'unit')
print("=========== RML Read ==========")
if apply_mapping:

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RML_CL3.ipdata'))
    print("============ RML data grid ============")
    evaluate_ordering()
    group_elem_parent_term(57)
    print("======== Grouped by parent element #57 =========")
    grow_tree(57, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RML_CL4_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RML_CL2.ipdata'))
    print("============ RML data grid ============")
    evaluate_ordering()
    group_elem_parent_term(58)
    print("======== Grouped by parent element #58 =========")
    grow_tree(58, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RML_CL1_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RML_CL1.ipdata'))
    print("============ RML data grid ============")
    evaluate_ordering()
    group_elem_parent_term(59)
    print("======== Grouped by parent element #59 =========")
    grow_tree(59, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RML_CL2_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RML_CL4.ipdata'))
    print("============ RML data grid ============")
    evaluate_ordering()
    group_elem_parent_term(60)
    print("======== Grouped by parent element #60 =========")
    grow_tree(60, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RML_CL3_mapping')
    print("==== RML Grown ======")
else:
    make_data_grid(0, 5.4, False, 'data_grid', 'data_grid')
    # define_data_geometry('RML_CL4.txt')
    print("============ RML data grid ============")
    evaluate_ordering()
    group_elem_parent_term(24)
    print("======== Grouped by parent element #24 =========")
    grow_tree(24, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, True, 'RML_mapping')


# Growing into RUL
define_node_geometry_2d(os.path.join(subject_path, 'Lung', 'SurfaceFEMesh', 'Right_fitted.ipnode'))
# define_node_geometry_2d('/hpc/bsha219/lung/Data/CTEPH/NewCase/FRC/Lung/SurfaceFEMesh/Right_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/RUL_surface.ipelem', 'unit')
print("=========== RUL mesh read ==========")
if apply_mapping:
    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RUL_CL1.ipdata'))
    print("============ RUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(35)
    print("======== Grouped by parent element #35 =========")
    grow_tree(35, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RUL_CL6_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RUL_CL5.ipdata'))
    print("============ RUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(36)
    print("======== Grouped by parent element #36 =========")
    grow_tree(36, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RUL_CL2_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RUL_CL6.ipdata'))
    print("============ RUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(51)
    print("======== Grouped by parent element #51 =========")
    grow_tree(51, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RUL_CL1_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RUL_CL2.ipdata'))
    print("============ RUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(52)
    print("======== Grouped by parent element #52 =========")
    grow_tree(52, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RUL_CL4_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RUL_CL3.ipdata'))
    print("============ RUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(49)
    print("======== Grouped by parent element #49 =========")
    grow_tree(49, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RUL_CL3_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'RUL_CL4.ipdata'))
    print("============ RUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(50)
    print("======== Grouped by parent element #50 =========")
    grow_tree(50, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'RUL_CL5_mapping')
    print("==== RUL Grown ======")
else:
    make_data_grid(0, 4.9, False, 'data_grid', 'data_grid')
    print("============ RUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(15)
    print("======== Grouped by parent element #15 =========")
    grow_tree(15, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, True, 'RUL_mapping')
    print("==== RUL Grown ======")


# Growing into LUL
define_node_geometry_2d(os.path.join(subject_path, 'Lung', 'SurfaceFEMesh', 'Left_fitted.ipnode'))
# define_node_geometry_2d('/hpc/bsha219/lung/Data/CTEPH/NewCase/FRC/Lung/SurfaceFEMesh/Left_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/LUL_surface.ipelem', 'unit')
print("=========== LUL Read ==========")
if apply_mapping:

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LUL_CL2.ipdata'))
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(27)
    print("======== Grouped by parent element #27 =========")
    grow_tree(27, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LUL_CL8_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LUL_CL6.ipdata'))
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(28)
    print("======== Grouped by parent element #28 =========")
    grow_tree(28, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LUL_CL7_mapping')
    print("==== LUL Grown ======")

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LUL_CL1.ipdata'))
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(29)
    print("======== Grouped by parent element #29 =========")
    grow_tree(29, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LUL_CL3_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LUL_CL3.ipdata'))
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(30)
    print("======== Grouped by parent element #30 =========")
    grow_tree(30, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LUL_CL6_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LUL_CL8.ipdata'))
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(41)
    print("======== Grouped by parent element #41 =========")
    grow_tree(41, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LUL_CL2_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LUL_CL7.ipdata'))
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(42)
    print("======== Grouped by parent element #42 =========")
    grow_tree(42, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, True, 'LUL_CL4_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LUL_CL4.ipdata'))
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(43)
    print("======== Grouped by parent element #43 =========")
    grow_tree(43, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LUL_CL1_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LUL_CL5.ipdata'))
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(44)
    print("======== Grouped by parent element #44 =========")
    grow_tree(44, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LUL_CL5_mapping')

else:
    make_data_grid(0, 5.3, False, 'data_grid', 'data_grid')
    print("============ LUL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(11)
    print("======== Grouped by parent element #11 =========")
    grow_tree(11, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, True, 'LUL_mapping')
    print("==== LUL Grown ======")

# Growing into LLL
define_node_geometry_2d(os.path.join(subject_path, 'Lung', 'SurfaceFEMesh', 'Left_fitted.ipnode'))
# define_node_geometry_2d('/hpc/bsha219/lung/Data/CTEPH/NewCase/FRC/Lung/SurfaceFEMesh/Left_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/LLL_surface.ipelem', 'unit')
print("=========== LLL Read ==========")
if apply_mapping:

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LLL_CL2.ipdata'))
    print("============ LLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(46)
    print("======== Grouped by parent element #46 =========")
    grow_tree(46, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LLL_CL3_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LLL_CL1.ipdata'))
    print("============ LLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(45)
    print("======== Grouped by parent element #45 =========")
    grow_tree(45, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LLL_CL2_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LLL_CL3.ipdata'))
    print("============ LLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(61)
    print("======== Grouped by parent element #61 =========")
    grow_tree(61, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LLL_CL1_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LLL_CL4.ipdata'))
    print("============ LLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(62)
    print("======== Grouped by parent element #62 =========")
    grow_tree(62, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LLL_CL4_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LLL_CL5.ipdata'))
    print("============ LLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(63)
    print("======== Grouped by parent element #63 =========")
    grow_tree(63, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LLL_CL6_mapping')

    define_data_geometry(os.path.join(subject_path, 'Intensity_mapping', 'LLL_CL6.ipdata'))
    print("============ LLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(64)
    print("======== Grouped by parent element #64 =========")
    grow_tree(64, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'LLL_CL5_mapping')
else:
    make_data_grid(0, 4.3, False, 'data_grid', 'data_grid')
    print("============ LLL data grid ============")
    evaluate_ordering()
    group_elem_parent_term(20)
    print("======== Grouped by parent element #12 =========")
    grow_tree(20, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, True, 'LLL_mapping')
    print("==== LLL Grown ======")


export_node_geometry('/hpc/bsha219/lung/Data/CTEPH/'+ subject + '/FRC/Vessel/' + subject + '_Artery_Full.exnode', 'MAC')
export_1d_elem_geometry('/hpc/bsha219/lung/Data/CTEPH/'+ subject + '/FRC/Vessel/' + subject + '_Artery_Full.exelem', 'MAC')
