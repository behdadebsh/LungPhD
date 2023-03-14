from aether.diagnostics import set_diagnostics_on
from aether.indices import define_problem_type
from aether.geometry import define_node_geometry_2d, define_elem_geometry_2d, make_data_grid, evaluate_ordering, \
    define_data_geometry, define_node_geometry, define_1d_elements, group_elem_parent_term
from aether.exports import export_node_geometry_2d, export_elem_geometry_2d, export_data_geometry, export_node_geometry, \
    export_1d_elem_geometry
from aether.growtree import grow_tree
from time import sleep

set_diagnostics_on(False)
define_problem_type('grow_tree')

define_node_geometry('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FTC/Vessel/CTEPH10_UpperArtery.ipnode')
define_1d_elements('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/HNA_A03_UpperArtery.ipelem')

# Growing into RLL
define_node_geometry_2d('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lung/SurfaceFEMesh/Right_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/lung/Data/UCSDLungs/HNA_A03/Supine/Vessel/RLL_surface.ipelem', 'unit')
make_data_grid(0, 6.0, False, 'test', 'test')
# export_data_geometry('wow','test',0)
# sleep(5)
evaluate_ordering()
group_elem_parent_term(23)
grow_tree(23, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'test')
# sleep(5)

# Growing into RUL
define_node_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Lung/SurfaceFEMesh/Right_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/RUL_surface.ipelem', 'unit')
make_data_grid(0, 6.0, False, 'test', 'test')
evaluate_ordering()
group_elem_parent_term(15)
grow_tree(15, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'test')

# Growing into RML
define_node_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Lung/SurfaceFEMesh/Right_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/RML_surface.ipelem', 'unit')
make_data_grid(0, 6.0, False, 'test', 'test')
evaluate_ordering()
group_elem_parent_term(24)
grow_tree(24, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'test')


# Growing into LUL
define_node_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Lung/SurfaceFEMesh/Left_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/LUL_surface.ipelem', 'unit')
make_data_grid(0, 6.0, False, 'test', 'test')
evaluate_ordering()
group_elem_parent_term(11)
grow_tree(11, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'test')

# Growing into LLL
define_node_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Lung/SurfaceFEMesh/Left_fitted.ipnode')
define_elem_geometry_2d('/hpc/bsha219/Hari/Data/UCSDLungs/HNA_A03/Supine/Vessel/LLL_surface.ipelem', 'unit')
make_data_grid(0, 6.0, False, 'test', 'test')
evaluate_ordering()
group_elem_parent_term(12)
grow_tree(12, 1, 60.0, 20.0, 0.4, 1.5, 1.5, 180.0, False, 'test')


export_node_geometry('Grown_Full.exnode', 'MAC')
export_1d_elem_geometry('Grown_Full.exelem', 'MAC')
