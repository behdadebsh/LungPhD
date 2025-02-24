import math
from collections import defaultdict
import statistics
import numpy as np
from numpy.linalg import norm
import os
import sys
import pymesh
import pyvista as pv

#############################################################################################################

#############################################################################################################

def get_connected_points(surf):
    # Extract the unique vertex IDs from the cells
    list_pts = set()
    for cell_id in range(surf.n_cells):
        list_pts.update(surf.cell_point_ids(cell_id))

    # Get the points associated with the unique vertex IDs
    points_surf = [surf.points[pnt] for pnt in list_pts]

    return points_surf

#############################################################################################################

# find the closest points in two point groups and return the average of their coordinates
def get_closest_between(gr1_points, gr2_points):
    closest_point = np.zeros((len(gr1_points)), dtype=int)
    shortest = np.full((len(gr1_points)), 1.0e10)
    
    for i, p1 in enumerate(gr1_points):
        for j, p2 in enumerate(gr2_points):
            distance = np.linalg.norm(p1 - p2)
            if distance < shortest[i]:
                shortest[i] = distance
                closest_point[i] = j
                
    ind_min = np.argmin(shortest)
    j_closest = closest_point[ind_min]
    
    between_point = (gr1_points[ind_min] + gr2_points[j_closest]) / 2.0
    
    return between_point

#############################################################################################################

def furthest_coord_from_point(point, coord):
    furthest_data = max((p for p in coord if p[1] < point[1]), key=lambda p: np.linalg.norm(point - p), default=None)
    return furthest_data

#############################################################################################################

def swap_smallest_from_a_to_b(surf, cell_ids1, cell_ids2):

    ids_swap = get_smallest_unconnected(surf, cell_ids1)
    for ids in ids_swap:
        cell_ids2.append(ids)
        cell_ids1.remove(ids)

    return cell_ids1, cell_ids2

#############################################################################################################

def swap_smallest_between(surf, cell_ids1, cell_ids2):

    ids_swap = get_smallest_unconnected(surf, cell_ids1)
    for ids in ids_swap:
        cell_ids2.append(ids)
        cell_ids1.remove(ids)

    ids_swap = get_smallest_unconnected(surf, cell_ids2)
    for ids in ids_swap:
        cell_ids1.append(ids)
        cell_ids2.remove(ids)

    return cell_ids1, cell_ids2

#############################################################################################################

def get_line_key(p0, p1):
    return str(min(p0, p1)) + '_' + str(max(p0, p1))

#############################################################################################################

def get_smallest_unconnected(surf, cell_ids):
    
    # given a surface with a selection of cells (defined by list of cell_ids), find all connected groups
    # and return the cells that are NOT in the largest connected group

    list_of_cells = []
    
    line_adj = defaultdict(list)
    for cell in cell_ids:
        ids = surf.cell_point_ids(cell)
        line_key = get_line_key(ids[0], ids[1])
        line_adj[line_key].append(cell)
        line_key = get_line_key(ids[0], ids[2])
        line_adj[line_key].append(cell)
        line_key = get_line_key(ids[2], ids[1])
        line_adj[line_key].append(cell)

    cell_adj = defaultdict(list)
    for key in line_adj:
        if len(line_adj[key]) > 1:
            cell1 = line_adj[key][0]
            cell2 = line_adj[key][1]
            cell_adj[cell1].append(cell2)
            cell_adj[cell2].append(cell1)

    len_groups = []
    ids_check = set(cell_ids)
    ids_check = list(ids_check)
    connected = defaultdict(list)
    nconnect = 0
    while len(ids_check) > 0:
        cell = ids_check[0]
        connected[nconnect].append(cell)
        layer = [cell]
        while len(layer) > 0: # while still some to check
            temp_layer = []  # record the next layer of cells
            for cell in layer:  # check the next layer of cells
                ids_check.remove(cell)
                for adj in cell_adj[cell]:
                    if adj not in connected[nconnect]:
                        connected[nconnect].append(adj)
                        temp_layer.append(adj)
            layer = list(temp_layer)
        len_groups.append(len(connected[nconnect]))
        nconnect +=1

    if nconnect > 1:
        ind_max = len_groups.index(max(len_groups))
        for i in range(nconnect):
            if i != ind_max:
                list_of_cells = list_of_cells + connected[i]

    #print ("number of connected groups=", nconnect)
    #print ("lengths", len_groups)
    #print ("swap", len(list_of_cells))
    return list_of_cells
            
#############################################################################################################

def output_mesh(name, vertices):
    write_data(name, vertices)
    export_data(name, vertices)
    
    
#############################################################################################################
    
def write_data(name, vertices):
    opf = open(name + '.ipdata', 'w+')
    opf.write("vertices\n")

    for nv in range(len(vertices)):
        opf.write("   %12d   %10.3f  %10.3f  %10.3f   1.0   1.0   1.0\n" % (nv+1, vertices[nv][0], vertices[nv][1], vertices[nv][2]))

    opf.close()
    

#############################################################################################################

def export_data(name, vertices):
    opf = open(name + '.exdata', 'w+')
    opf.write(" Group name: for_fitting\n")
    opf.write(" #Fields=1\n")
    opf.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
    opf.write("   x.  Value index= 1, #Derivatives=0\n")
    opf.write("   y.  Value index= 2, #Derivatives=0\n")
    opf.write("   z.  Value index= 3, #Derivatives=0\n")

    for nv in range(len(vertices)):
        opf.write(" Node: %12d\n" % (nv+1))
        opf.write("   %10.3f  %10.3f  %10.3f\n" % (vertices[nv][0], vertices[nv][1], vertices[nv][2]))

    opf.close()
    
#############################################################################################################

# assemble a list of all triangle lines
def get_list_of_lines_ids(surf, ids_surf):
    list_lines = []
    for cell in ids_surf:
        pts = surf.cell_point_ids(cell)
        line_key = get_line_key(pts[0], pts[1])
        list_lines.append(line_key)
        line_key = get_line_key(pts[0], pts[2])
        list_lines.append(line_key)
        line_key = get_line_key(pts[2], pts[1])
        list_lines.append(line_key)

    return list_lines


#############################################################################################################

# assemble a list of all triangle points

def get_list_of_points_ids(surf, ids_surf):

    list_points = []
    for cell in ids_surf:
        pts = surf.cell_point_ids(cell)
        list_points.append(pts[0])
        list_points.append(pts[1])
        list_points.append(pts[2])

    list_points = list(set(list_points))
    return list_points


#############################################################################################################

# Get the distance to the closest point in gr_points from the given point
def get_shortest_dist_to_point(point, gr_points):
    distances = [np.linalg.norm(np.array(point) - np.array(p)) for p in gr_points]
    shortest = min(distance)
    return shortest

#############################################################################################################

# Get the index of the closest point in gr_points to the given point
def get_closest_point_to_group(point, gr_points):
    distances = [np.linalg.norm(np.array(point) - np.array(p)) for p in gr_points]
    j = np.argmin(distances)
    return j

#############################################################################################################

def plane_from_3pts(p0, p1, p2):

    vect_v = np.array(p0 - p2)
    vect_u = np.array(p1 - p2)
    normal = np.cross(vect_u, vect_v)
    normal = normal/np.linalg.norm(normal)
    d = -np.dot(p0, normal)

    return normal, d

#############################################################################################################

def should_swap(cell, y, z, ygrt, zgrt, y_lim, z_lim):
    if ygrt and zgrt:
        return y > y_lim and z > z_lim
    elif zgrt:
        return y < y_lim and z > z_lim
    else:
        return y > y_lim and z < z_lim

#############################################################################################################

def swap_cells_using_plane(surf, cell_centers, ids1, ids2, normal, d, y_lim, z_lim, ygrt, zgrt):

    swap_ids = []
    for cell in ids1:
        y, z = cell_centers.points[cell, 1], cell_centers.points[cell, 2]
        if should_swap(cell, y, z, ygrt, zgrt, y_lim, z_lim):
            dist = np.dot(cell_centers.points[cell], normal) + d
            if (dist > 0.0 and ygrt and zgrt) or (dist < 0.0 and not ygrt and not zgrt):
                swap_ids.append(cell)

    for cell in swap_ids:
        ids2.append(cell)
        ids1.remove(cell)

    swap_ids = []
    for cell in ids2:
        y, z = cell_centers.points[cell, 1], cell_centers.points[cell, 2]
        if should_swap(cell, y, z, ygrt, zgrt, y_lim, z_lim):
            dist = np.dot(cell_centers.points[cell], normal) + d
            if (dist < 0.0 and ygrt and zgrt) or (dist > 0.0 and not ygrt and not zgrt):
                swap_ids.append(cell)

    for cell in swap_ids:
        ids1.append(cell)
        ids2.remove(cell)

    return ids1, ids2

#############################################################################################################

# swap cell ids from list a to list b based on the distance from a user-defined clipping plane

def swap_cells_from_a_to_b(cell_centers, ids1, ids2, normal, d, y_min, y_max, z_min, invert):
    # cell_centers == coordinates of the centres of surface cells
    # ids1 == list of cell ids for group to swap FROM
    # ids2 == list of cell ids for group to swap TO
    # normal, d == defines the clipping plane
    # y_min, y_max, z_min == defines the range over which swapping occurs
    # invert == defines which 'side' of plane the swapping applies to
    
    swap_ids = [] # accumulates list of cell ids to swap
    for cell in ids1:
        y, z = cell_centers.points[cell, 1], cell_centers.points[cell, 2]
        if y < y_max and y > y_min and z > z_min:
            dist = np.dot(cell_centers.points[cell], normal) + d # distance to plane
            if (not invert and dist < 0.0) or (invert and dist > 0.0):
                swap_ids.append(cell)

    for cell in swap_ids:
        ids2.append(cell)
        ids1.remove(cell)

    return ids1, ids2

#############################################################################################################

# identify the points on the apical lobe (left or right) that are closest to the 'inflexion' and 'anterior'
# identified for the purpose of using clipping planes to reallocate apical 'knotch' points from the lateral
# to the medial surface, and to smooth the apical region.

def get_anterior_inflexion(COM_UL, init_surface, ids, pct_ant, pct_inf):
    # COM_UL == centre of mass of the upper lobe surfaces (which this is written for)
    # init_surface == entire upper lobe surface
    # ids == list of cell ids to get just one annoted 'side' of init_surface (medial or lateral)
    # pct_ant == location of the anterior point along z-axis of lobe
    # pct_inf == location of the inflexion point along z-axis of lobe

    # extract surface cells and edge for the input surface list (ids)
    sub_surface = init_surface.extract_cells(ids).extract_largest()
    sub_edge = edge_points(sub_surface)

    # calculate max, min coordinate bounds for the surface cells
    bounds = sub_surface.bounds
    y_threshold = 0.3*bounds[2] + 0.7*bounds[3] # hard-coded; 'pushes' targets to anterior
    target_z1 = (1.0 - pct_inf)*bounds[4] + pct_inf*bounds[5] # target z-coord for inflexion point
    target_z2 = (1.0 - pct_ant)*bounds[4] + pct_ant*bounds[5] # target z-coord for anterior point

    # Filter points with y < y_threshold
    filtered_points = sub_edge[sub_edge[:, 1] < y_threshold] # remove posterior points

    # Find the index of the point with the z coordinate closest to the target_z
    index_closest_z1 = np.argmin(np.abs(filtered_points[:, 2] - target_z1)) # for inflexion
    index_closest_z2 = np.argmin(np.abs(filtered_points[:, 2] - target_z2)) # for anterior

    # Get the closest point's coordinates
    inflex_point = filtered_points[index_closest_z1]
    anteri_point = filtered_points[index_closest_z2]

    return anteri_point, inflex_point
    
#############################################################################################################

# identify 'necks' as narrow regions on a surface

def find_neck(surf, ids_surf0, ids_surf1, ids_surf2):

    list_points0 = get_list_of_points_ids(surf, ids_surf0)
    list_points1 = get_list_of_points_ids(surf, ids_surf1)
    list_points2 = get_list_of_points_ids(surf, ids_surf2)
    swap_cells = []
    for cell in ids_surf0: # the surface we are checking for a neck
        pts = surf.cell_point_ids(cell)
        count1 = list_points1.count(pts[0]) + list_points2.count(pts[1])
        count2 = list_points1.count(pts[0]) + list_points2.count(pts[2])
        count3 = list_points1.count(pts[2]) + list_points2.count(pts[1])
        count4 = list_points2.count(pts[0]) + list_points1.count(pts[1])
        count5 = list_points2.count(pts[0]) + list_points1.count(pts[2])
        count6 = list_points2.count(pts[2]) + list_points1.count(pts[1])
        if count1 == 2 or count2 == 2 or count3 == 2 or count3 == 2 or count4 == 2 or count5 == 2 or count6 == 2:
            swap_cells.append(cell)
    
    return swap_cells

#############################################################################################################

# smooth the edge between two given surfaces (defined by id lists ids_surf1, ids_surf2) by swapping cells
# that have two or more edges in the 'other' surface into that surface

def smooth_edge(surf, ids_surf1, ids_surf2):

    count = 0
    n = 1
    m = 1
    while (n > 0 or m > 0) and count < 5:
        n, ids_swap1, ids_swap2 = check_for_swaps(surf, ids_surf1, ids_surf2)
        m, ids_swap2, ids_swap1 = check_for_swaps(surf, ids_swap2, ids_swap1)

    return ids_swap1, ids_swap2
    
#############################################################################################################

# check all triangles to identify ones that do not belong in the group. get the line ids for all lines in the
# first surface. Then count the number of times that lines for each cell from the second surface appear in the
# first (list of lines). If more than once, then swap surfaces.

def check_for_swaps(surf, ids_surf1, ids_surf2):

    list_lines = get_list_of_lines_ids(surf, ids_surf1)
    swap_cells = []
    for cell in ids_surf2:
        pts = surf.cell_point_ids(cell)
        line1 = get_line_key(pts[0], pts[1])
        line2 = get_line_key(pts[0], pts[2])
        line3 = get_line_key(pts[2], pts[1])
        count_lines = list_lines.count(line1) + list_lines.count(line2) + list_lines.count(line3)
        if count_lines > 1:
            swap_cells.append(cell)

    for cell in swap_cells:
        ids_surf1.append(cell)
        ind_surf2 = ids_surf2.index(cell)
        del ids_surf2[ind_surf2]

    return len(swap_cells), ids_surf1, ids_surf2
    
    
#############################################################################################################

# check all triangles to identify ones that do not belong in the group based on their points. 

def check_for_swaps_points(surf, ids_surf1, ids_surf2):

    list_points = get_list_of_points_ids(surf, ids_surf1) # point ids in group ids_surf1
    swap_cells = []
    for cell in ids_surf2:
        pts = surf.cell_point_ids(cell)
        count_points = list_points.count(pts[0]) + list_points.count(pts[1]) + list_points.count(pts[2])
        if count_points == 3:
            swap_cells.append(cell)

    for cell in swap_cells:
        ids_surf1.append(cell)
        ind_surf2 = ids_surf2.index(cell)
        del ids_surf2[ind_surf2]

    return ids_surf1, ids_surf2
    
    
#############################################################################################################

def fix_mesh(mesh, detail):
    bbox_min, bbox_max = mesh.bbox
    diag_len = norm(bbox_max - bbox_min)
    if detail == "normal":
        target_len = diag_len * 5e-3
    elif detail == "high":
        target_len = diag_len * 2.5e-3
    elif detail == "low":
        target_len = diag_len * 1e-2
    elif detail == "extra_low":
        target_len = diag_len * 0.02
    #print("Target resolution: {} mm".format(target_len))

    count = 0
    # several lines commented out because while they do improve the 'mesh'
    # they are not needed for the purposes of reducing the number of
    # data points
    #mesh, __ = pymesh.remove_degenerated_triangles(mesh, 100)
    #mesh, __ = pymesh.split_long_edges(mesh, target_len)
    num_vertices = mesh.num_vertices
    while True:
        mesh, __ = pymesh.collapse_short_edges(mesh, 1e-6)
        mesh, __ = pymesh.collapse_short_edges(mesh, target_len,
                                               preserve_feature=True)
        mesh, __ = pymesh.remove_obtuse_triangles(mesh, 150.0, 100)
        if mesh.num_vertices == num_vertices:
            break

        num_vertices = mesh.num_vertices
        #print("#v: {}".format(num_vertices))
        count += 1
        if count > 10: break

    mesh = pymesh.resolve_self_intersection(mesh)
    mesh, __ = pymesh.remove_duplicated_faces(mesh)
    mesh = pymesh.compute_outer_hull(mesh)
    mesh, __ = pymesh.remove_duplicated_faces(mesh)
    mesh, __ = pymesh.remove_obtuse_triangles(mesh, 179.0, 5)
    mesh, __ = pymesh.remove_isolated_vertices(mesh)

    return mesh, target_len


#############################################################################################################

def is_smallest(value, lst):
    return value == min(lst)

#############################################################################################################

def angle_to_xyz(dirn, normal):
    directions = { "x+" : np.array([1, 0, 0]), "x-" : np.array([-1, 0, 0]), "y+" : np.array([0, 1, 0]),
                   "y-" : np.array([0, -1, 0]), "z+" : np.array([0, 0, 1]), "z-" : np.array([0, 0, -1]) }
    value = math.acos(np.dot(directions[dirn], normal))
    return value

#############################################################################################################

def ray_through_one_cell(mesh, point0, point1, cell):
    value = False
    point, index = mesh.ray_trace(point0, point1)
    if len(index) == 0 or (len(index) == 1 and index[0] == cell):
        value = True
    return value

#############################################################################################################

def closest_edge_points(edge1, edge2):

    edge1_points = edge1.points
    edge2_points = edge2.points
    min_distance = 1.0e+10
    
    for i, point in enumerate(edge1_points):
        index = edge2.find_closest_point(point)
        distance = math.dist(point, edge2_points[index])
        if distance < min_distance:
            index1 = i
            index2 = index
            min_distance = distance
    return edge1_points[index1], edge2_points[index2]
   
#############################################################################################################

def edge_points(mesh):
    edges = mesh.extract_feature_edges(boundary_edges=True, non_manifold_edges=True, feature_edges=False, manifold_edges=False)
    values = edges.points
    return values

#############################################################################################################

def annotate_edge(points1, points2, dist_lim):
    common_pts = np.empty((0, 3))
    for i, point1 in enumerate(points1):
        for j, point2 in enumerate(points2):
            distance = math.dist(point1, point2)
            if distance < dist_lim:
                common_pts = np.vstack((common_pts, point1))
    return common_pts

#############################################################################################################

def right_surface_data(output_directory, phenotype, dist_lim, sharp):

    RLL_surf = pv.read(os.path.join(output_directory, 'RLL_surf.stl'))
    RML_surf = pv.read(os.path.join(output_directory, 'RML_surf.stl'))
    RUL_surf = pv.read(os.path.join(output_directory, 'RUL_surf.stl'))

    #---- separate the fissures from the lobes
    
    # identify and remove the RH fissure from the RUL surface
    _ = RUL_surf.compute_implicit_distance(RLL_surf, inplace=True)
    RU_L_fiss_surface = RUL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)

    if sharp:
        RUL_rem = RUL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False)
        RUL_rem.compute_implicit_distance(RML_surf, inplace=True)
        RUL_rem['implicit_distance'] = np.abs(RUL_rem['implicit_distance'])
        RU_M_fiss_surface = RUL_rem.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)
        RUL_rem = RUL_rem.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False)
    else:
        _ = RUL_surf.compute_implicit_distance(RML_surf, inplace=True)
        RU_M_fiss_surface = RUL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True) #, inplace=True)
        R_horiz_fiss = RU_L_fiss_surface.merge(RU_M_fiss_surface)
        points = pv.wrap(R_horiz_fiss.points)
        RHF_merge = points.reconstruct_surface()
        _ = RUL_surf.compute_implicit_distance(RHF_merge, inplace=True)
        max_dist = max(RUL_surf.point_data['implicit_distance'])
        min_dist = min(RUL_surf.point_data['implicit_distance'])
        if abs(max_dist) > abs(min_dist):
            RUL_rem = RUL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False) #, inplace=True)
        else:
            RUL_rem = RUL_surf.clip_scalar(scalars='implicit_distance', value = -dist_lim, invert=True) #, inplace=True)

    RUL_rem = RUL_rem.extract_largest()
    RUL_edges = RUL_rem.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, feature_edges=False, manifold_edges=False)
    R_horiz_fiss = RU_L_fiss_surface.merge(RU_M_fiss_surface)

    RLL_edges = RLL_surf.extract_feature_edges(feature_angle=60.0, feature_edges=True)
    # identify and remove the RO fissure from the RLL surface
    RLL_surf.compute_implicit_distance(RUL_surf, inplace=True)
    RL_U_fiss_surface = RLL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)
    
    if sharp:
        RLL_rem = RLL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False)
        RLL_rem.compute_implicit_distance(RML_surf, inplace=True)
        RLL_rem['implicit_distance'] = np.abs(RLL_rem['implicit_distance'])
        RL_M_fiss_surface_max = RLL_rem.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)
        RL_M_fiss_surface = RL_M_fiss_surface_max.clip_scalar(scalars='implicit_distance', value = -dist_lim, invert=False)
        RLL_rem = RLL_rem.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False)
    else:
        _ = RLL_surf.compute_implicit_distance(RML_surf, inplace=True)
        RL_M_fiss_surface = RLL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)
        R_oblique_fiss = RL_U_fiss_surface.merge(RL_M_fiss_surface)
        points = pv.wrap(R_oblique_fiss.points)
        ROF_merge = points.reconstruct_surface()
        _ = RLL_surf.compute_implicit_distance(ROF_merge, inplace=True)
        RLL_rem = RLL_surf.clip_scalar(scalars='implicit_distance', value = -dist_lim, invert=True)

    RLL_rem = RLL_rem.extract_largest()
    R_oblique_fiss = RL_U_fiss_surface.merge(RL_M_fiss_surface)
    
    # identify and remove the RO and RH fissures from the RML surface
    _ = RML_surf.compute_implicit_distance(RUL_surf, inplace=True)
    RM_horiz_fiss = RML_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)
    RML_rem = RML_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False)
    _ = RML_rem.compute_implicit_distance(RLL_surf, inplace=True)
    RM_oblique_fiss = RML_rem.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)
    RML_rem = RML_rem.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False)

    RML_rem = RML_rem.extract_largest().extract_surface()
    mean_surface_area = RML_rem.area/RML_rem.n_cells
    average_length = math.sqrt((4*mean_surface_area)/math.sqrt(3))

    COM_RHM = RM_horiz_fiss.center_of_mass()
    COM_ROM = RM_oblique_fiss.center_of_mass()
    COM_RL_M_fiss = RL_M_fiss_surface.center_of_mass()

    # ---- RIGHT LOWER LOBE ----
    bounds_coords = RLL_rem.bounds
    bounds_centre = RLL_rem.center
    base_centre = np.array([bounds_centre[0], bounds_centre[1], bounds_coords[4]])

    av_fiss_normal = RL_M_fiss_surface.compute_normals(cell_normals=True)['Normals'].mean(axis=0)
    av_fiss_normal /= np.linalg.norm(av_fiss_normal)

    cell_normals = RLL_rem.compute_normals(cell_normals=True)['Normals']
    cell_centers = RLL_rem.cell_centers()
    ids_bas, ids_lat, ids_med, ids_fiss = [], [], [], []

    for cell in range(RLL_rem.n_cells):
        cell_normal = cell_normals[cell]
        cell_centre = cell_centers.points[cell]
        angles = [angle_to_xyz("x-", cell_normal), angle_to_xyz("x+", cell_normal),  
                  angle_to_xyz("y+", cell_normal), angle_to_xyz("z-", cell_normal)]
        fiss_angle = math.acos(np.dot(cell_normal, av_fiss_normal))
        if cell_centre[2] < COM_RL_M_fiss[2] and ray_through_one_cell(RLL_rem, cell_centre, base_centre, cell):
            ids_bas.append(cell)
        elif cell_centre[0] > COM_RL_M_fiss[0] and is_smallest(angles[3], angles):
            ids_bas.append(cell)
        elif cell_centre[0] <= COM_RL_M_fiss[0] and is_smallest(angles[3], angles) and angles[3] < fiss_angle:
            ids_bas.append(cell)
        elif is_smallest(angles[0], angles[0:1]):
            ids_lat.append(cell)
        else:
            ids_med.append(cell)

    # swap isolated cells between the two groups. this smoothes up the edges a bit
    ids_bas, ids_med = smooth_edge(RLL_rem, ids_bas, ids_med) # only swap into the base
    ids_lat, ids_bas = smooth_edge(RLL_rem, ids_lat, ids_bas) # only swap out to the lateral
    ids_med, ids_bas = smooth_edge(RLL_rem, ids_med, ids_bas) # only swap out to the medial

    # RLL base is complete. Extract the largest connected component as the RLL base
    RLL_bas = RLL_rem.extract_cells(ids_bas).extract_largest()

    # find the closest points on the RLL base edge and R oblique fissure edge
    RO_fiss_edges = R_oblique_fiss.extract_feature_edges(boundary_edges=True) 
    RLL_base_edges = RLL_bas.extract_feature_edges(boundary_edges=True)
    closest_RLF, closest_RLB = closest_edge_points(RO_fiss_edges, RLL_base_edges)
    
    # CORRECT ALLOCATION OF POINTS TOWARDS RLL ANTERIOR
    # get the cells that fall on each side of a dividing plane, based
    # on location of fissure anterior base point, COM of the base and COM of the RO fissure
    COM_RL_base = RLL_bas.center_of_mass()
    COM_RO_fiss = R_oblique_fiss.center_of_mass()
    normal, d = plane_from_3pts(closest_RLF, COM_RL_base, COM_RO_fiss)
    y_lim = 0.2*bounds_coords[2] + 0.8*bounds_coords[3] 
    z_lim = closest_RLB[2] - 30.0
    ids_lat, ids_med  = swap_cells_using_plane(RLL_rem, cell_centers, ids_lat, ids_med,
                                              normal, d, y_lim, z_lim, False, True)

    # CORRECT ALLOCATION OF POINTS TOWARDS RLL POSTERIOR
    post_base_RLL = RLL_bas.points[np.argmax(RLL_bas.points[:, 1])] # location on base of max y
    post_apex_RLL = R_horiz_fiss.points[np.argmax(R_horiz_fiss.points[:, 1])] # location on RH fissure of max y

    # get the cells that fall on each side of a dividing plane, based
    # on location of apical, and posterior RLL apical points
    COM_RLL = RLL_rem.center_of_mass()
    normal, d = plane_from_3pts(post_apex_RLL, post_base_RLL, COM_RLL)
    y_lim = 0.8*bounds_coords[2] + 0.2*bounds_coords[3] 
    z_lim = post_base_RLL[2] - 30
    ids_lat, ids_med  = swap_cells_using_plane(RLL_rem, cell_centers, ids_lat, ids_med,
                                              normal, d, y_lim, z_lim, True, True)
    # smooth out a bit
    ids_lat, ids_med = smooth_edge(RLL_rem, ids_lat, ids_med) # to lateral
    #ids_lat, ids_med = check_for_swaps(RLL_rem, ids_lat, ids_med) # to lateral
    #ids_med, ids_lat = check_for_swaps(RLL_rem, ids_med, ids_lat) # to medial
    ids_lat, ids_med = swap_smallest_from_a_to_b(RLL_rem, ids_lat, ids_med) # from lateral to medial
    ids_bas, ids_med = swap_smallest_from_a_to_b(RLL_rem, ids_bas, ids_med) # from basal to medial

    RLL_med = RLL_rem.extract_cells(ids_med).extract_largest()
    RLL_lat = RLL_rem.extract_cells(ids_lat).extract_largest()

    RLL_bas_edge = edge_points(RLL_bas)
    RLL_med_edge = edge_points(RLL_med)
    RLL_lat_edge = edge_points(RLL_lat)
    RL_M_fis_edge = edge_points(RL_M_fiss_surface)
    RL_U_fis_edge = edge_points(RL_U_fiss_surface)

    edge_RL_base_lat = annotate_edge(RLL_bas_edge, RLL_lat_edge, average_length/2)
    edge_RL_base_med = annotate_edge(RLL_bas_edge, RLL_med_edge, average_length/2)
    edge_RL_fiss_base = annotate_edge(RLL_bas_edge, RL_M_fis_edge, average_length/2)
    edge_RL_posterior = annotate_edge(RLL_lat_edge, RLL_med_edge, average_length/2)
    edge_RL_U_fiss_lat = annotate_edge(RL_U_fis_edge, RLL_lat_edge, average_length/2)
    edge_RL_U_fiss_med = annotate_edge(RL_U_fis_edge, RLL_med_edge, average_length/2)
    edge_RL_M_fiss_lat = annotate_edge(RL_M_fis_edge, RLL_lat_edge, average_length/2)
    edge_RL_M_fiss_med = annotate_edge(RL_M_fis_edge, RLL_med_edge, average_length/2)

    # -------------------------------------------------------------------------- end RLL
    
    # ---- RIGHT UPPER LOBE ----
    apex_RUL = RUL_rem.points[np.argmax(RUL_rem.points, axis=0)[2]] 
    COM_RH_U = RU_M_fiss_surface.center_of_mass()
    ant_RH_U = furthest_coord_from_point(COM_RH_U, RU_M_fiss_surface.points)
    COM_RUL = RUL_rem.center_of_mass()
    normal, d = plane_from_3pts(apex_RUL, post_apex_RLL, COM_RUL)

    bounds_coords = RUL_rem.bounds
    y_lim = apex_RUL[1] - 0.25*(bounds_coords[3] - bounds_coords[2])
    z_lim = ant_RH_U[2] + 0.1*(bounds_coords[5] - bounds_coords[4])
    centre_RUL = np.array([bounds_coords[1], COM_RH_U[1], COM_RUL[2]])

    cell_centers = RUL_rem.cell_centers()
    ids_lat, ids_med = [], []

    # find the 'knee' on anterior, taken as the most apical point in the anterior region
    # (defined as y < y_lim) that will be classified as medial (using ray tracing)
    knee = np.array([-1.0e8,-1.0e8,-1.0e8])
    for cell in range(RUL_rem.n_cells):
        x, y, z = cell_centers.points[cell]
        if x > COM_RUL[0] and y < y_lim: #  using ray tracing to identify cells on RUL medial
            point, index = RUL_rem.ray_trace(cell_centers.points[cell], centre_RUL)
            if not index.any() or (len(index) == 1 and index[0] == cell):
                if z > knee[2]:
                    knee = (x, y, z)

    normal2, d2 = plane_from_3pts(apex_RUL, knee, COM_RUL)
    normal3, d3 = plane_from_3pts(ant_RH_U, ant_RH_U + np.array([50,0,100]), COM_RHM)
                    
    for cell in range(RUL_rem.n_cells):
        x,y,z = cell_centers.points[cell]
        if y > apex_RUL[1]:  # divide points using the plane that contains apex_RUL, post_apex_RLL, COM_RUL
            dist = np.dot(cell_centers.points[cell], normal) + d
            if dist < 0.0: 
                ids_lat.append(cell)
            else:
                ids_med.append(cell)
        elif x < COM_RH_U[0]:  # assume all cells on lateral side of COM are lateral surface
            ids_lat.append(cell)
        else:
            if y > knee[1]:  # for points that are more posterior than the knee point,
                if z > COM_RUL[2]:  # divide using the plane that contains apex_RUL, knee, COM_RUL
                    dist = np.dot(cell_centers.points[cell], normal2) + d2
                    if dist > 0.0: 
                        ids_lat.append(cell)
                    else:
                        ids_med.append(cell)
                else:
                    dist = np.dot(cell_centers.points[cell], normal3) + d3
                    if dist < 0.0:   # divide points using the plane that contains ant_RH, COM_RH fissure
                        ids_lat.append(cell)
                    else:
                        ids_med.append(cell)
            else:
                dist = np.dot(cell_centers.points[cell], normal3) + d3
                if z < z_lim and dist < 0.0:
                    ids_lat.append(cell)
                else:
                    #  using ray tracing to identify cells on RUL medial
                    z_centre = min(np.array([z, knee[2], 0.5*(knee[2]+COM_RUL[2])]))
                    centre_lungs = np.array([bounds_coords[1], knee[1], z_centre])
                    point, index = RUL_rem.ray_trace(cell_centers.points[cell], centre_lungs)
                    if len(index) == 0 or (len(index) == 1 and index[0] == cell):
                        ids_med.append(cell)
                    else:
                        ids_lat.append(cell)

    # find edge triangles and small unconnected groups that are mislabelled and swap to other group
    ids_med, ids_lat = swap_smallest_between(RUL_rem, ids_med, ids_lat)

    # clip apical 'knotch' using a clipping plane
    y_max = apex_RUL[1]
    frac_ant_pnt = 0.6  # location of the anterior point along z-axis of RUL
    frac_inf_pnt = 0.85  # location of the inflexion point along z-axis of RUL
    anterior_point, inflexion_point = get_anterior_inflexion(COM_RUL, RUL_rem, ids_lat, frac_ant_pnt, frac_inf_pnt)

    # clip between inflexion point and anterior
    normal, d = plane_from_3pts(COM_RUL, inflexion_point, anterior_point)
    y_min = anterior_point[1]
    z_min = anterior_point[2]
    ids_med, ids_lat = swap_cells_from_a_to_b(cell_centers, ids_med, ids_lat, normal, d, y_min, y_max, z_min, True)

    # clip between inflexion point and apical point
    normal, d = plane_from_3pts(COM_RUL, apex_RUL, inflexion_point)
    y_min = inflexion_point[1]
    z_min = inflexion_point[2]
    ids_lat, ids_med = swap_cells_from_a_to_b(cell_centers, ids_lat, ids_med, normal, d, y_min, y_max, z_min, False)

    ids_med, ids_lat = smooth_edge(RUL_rem, ids_med, ids_lat) # check for swaps into the medial group
    #ids_med, ids_lat = check_for_swaps(RUL_rem, ids_med, ids_lat) # check for swaps into the medial group
    #ids_lat, ids_med = check_for_swaps(RUL_rem, ids_lat, ids_med) # check for swaps into the lateral group
    
    RUL_med = RUL_rem.extract_cells(ids_med).extract_largest()
    RUL_lat = RUL_rem.extract_cells(ids_lat).extract_largest()

    RUL_med_edge = edge_points(RUL_med)
    RUL_lat_edge = edge_points(RUL_lat)
    RU_M_fis_edge = edge_points(RU_M_fiss_surface)
    RU_L_fis_edge = edge_points(RU_L_fiss_surface)

    edge_RU_ant_post = annotate_edge(RUL_lat_edge, RUL_med_edge, average_length/2)
    edge_RU_L_fiss_lat = annotate_edge(RU_L_fis_edge, RUL_lat_edge, average_length/2)
    edge_RU_L_fiss_med = annotate_edge(RU_L_fis_edge, RUL_med_edge, average_length/2)
    edge_RU_M_fiss_lat = annotate_edge(RU_M_fis_edge, RUL_lat_edge, average_length/2)
    edge_RU_M_fiss_med = annotate_edge(RU_M_fis_edge, RUL_med_edge, average_length/2)

    # -------------------------------------------------------------------------- end RUL
    
    # ---- RIGHT MIDDLE LOBE ----

    ## from here is the trial code
    COM_RML = RML_rem.center_of_mass()
    bounds_coords = RML_rem.bounds
    bounds_centre = RML_rem.center
    base_centre = np.array([bounds_centre[0], bounds_centre[1], bounds_coords[4]])
    ant_RH = get_closest_between(edge_RU_M_fiss_lat, edge_RU_M_fiss_med)
    hilum = get_closest_between(edge_RU_M_fiss_med, edge_RU_L_fiss_med)
    ant_base_lat_RLL = get_closest_between(edge_RL_base_lat, edge_RL_M_fiss_lat)
    ant_base_med_RLL = get_closest_between(edge_RL_base_med, edge_RL_M_fiss_med)
    cent_RL_fiss = np.mean(edge_RL_fiss_base, axis=0)
    normal, d = plane_from_3pts(ant_RH, cent_RL_fiss, COM_ROM)
    normal2, d2 = plane_from_3pts(ant_RH, hilum, bounds_centre)
    normal3, d3 = plane_from_3pts(ant_RH, ant_base_med_RLL, bounds_centre)
    #dist_ant_RH_to_hilum = math.dist(ant_RH, hilum)
    dist_ant_base_to_hilum = math.dist(ant_base_med_RLL, hilum)
    
    sRLL_bas = RLL_bas.extract_surface()
    sRUL_lat = RUL_lat.extract_surface()
    sRLL_med = RLL_med.extract_surface()
    sRUL_med = RUL_med.extract_surface()
    av_base_normal = sRLL_bas.compute_normals(cell_normals=True)['Normals'].mean(axis=0)
    av_base_normal /= np.linalg.norm(av_base_normal)
    av_med_normal = sRLL_med.compute_normals(cell_normals=True)['Normals'].mean(axis=0)
    av_med_normal /= np.linalg.norm(av_med_normal)
    base_centre = RLL_bas.center_of_mass()

    cell_centers = sRUL_lat.cell_centers()
    cell_normals = sRUL_lat.compute_normals(cell_normals=True)['Normals']
    av_lat_normal = [0, 0, 0]
    for cell in range(sRUL_lat.n_cells):
        x,y,z = cell_centers.points[cell]
        if y < bounds_centre[2]: #bounds_coords[3]:
            av_lat_normal = av_lat_normal + cell_normals[cell]
    av_lat_normal /= np.linalg.norm(av_lat_normal) # don't need to divide by N because normalising anyway
    lat_focus_pt = bounds_centre + av_lat_normal * 150.0
    base_pt = np.array(cent_RL_fiss) - np.array([0, 10, 50])
    #base_pt = np.array(cent_RL_fiss) - np.array([0, 20, 40])
    #base_pt = np.array(ant_base_lat_RLL) - np.array([-55, 25, 10])

    cell_centers = sRUL_med.cell_centers()
    cell_normals = sRUL_med.compute_normals(cell_normals=True)['Normals']
    av_med_normal = [0, 0, 0]
    for cell in range(sRUL_med.n_cells):
        x,y,z = cell_centers.points[cell]
        if y < bounds_centre[2]:
            av_med_normal = av_med_normal + cell_normals[cell]
    av_med_normal /= np.linalg.norm(av_med_normal) # don't need to divide by N because normalising anyway

    cell_normals = RML_rem.compute_normals(cell_normals=True)['Normals']
    cell_centers = RML_rem.cell_centers()
    ids_bas, ids_lat, ids_med, ids_med2 = [], [], [], []

    for cell in range(RML_rem.n_cells):
        cell_normal = cell_normals[cell]
        cell_centre = cell_centers.points[cell]
        x,y,z = cell_centers.points[cell]
        dist_to_hilum = math.dist(cell_centre, hilum)
        lat_angle = math.acos(np.dot(cell_normal, av_lat_normal))
        med_angle = math.acos(np.dot(cell_normal, av_med_normal))
        base_angle = math.acos(np.dot(cell_normal, av_base_normal))
        
        angles = [lat_angle, med_angle, base_angle]
        dist = np.dot(cell_centers.points[cell], normal) + d
        dist2 = np.dot(cell_centers.points[cell], normal2) + d2
        dist3 = np.dot(cell_centers.points[cell], normal3) + d3

        point, index = RML_rem.ray_trace(cell_centers.points[cell], base_pt)
        if not index.any() or (len(index) == 1 and index[0] == cell):
            single_base_ray = True
        else:
            single_base_ray = False

        bdrn = (base_centre - cell_centers.points[cell]) # use the angle to the base_centre
        bdrn = bdrn/np.linalg.norm(bdrn)
        b_angle = math.acos(np.dot(bdrn, cell_normal))
        if is_smallest(abs(lat_angle), [abs(lat_angle), abs(med_angle)]) and dist > -50.0:
            ids_lat.append(cell)
        elif dist3 < 0.0:
            ids_med.append(cell)
        elif single_base_ray:
            ids_bas.append(cell)
        else:
            ids_med2.append(cell) # magenta, good
    for cell in range(RML_rem.n_cells):
        if cell not in ids_lat and cell not in ids_med and cell not in ids_bas and cell not in ids_med2:
            print (cell)

    ids_lat, ids_med = swap_smallest_from_a_to_b(RML_rem, ids_lat, ids_med) # from lateral to medial
    ids_med2, ids_bas = swap_smallest_from_a_to_b(RML_rem, ids_med2, ids_bas) # from medial to basal
    ids_bas, ids_med2 = check_for_swaps_points(RML_rem, ids_bas, ids_med2) # from medial to base
    ids_bas, ids_lat = smooth_edge(RML_rem, ids_bas, ids_lat)
    ids_med += ids_med2 # append ids_med with ids_med2 (don't do earlier; previous swapping important)
    ids_med, ids_bas = swap_smallest_from_a_to_b(RML_rem, ids_med, ids_bas) # from medial to basal
    ids_bas, ids_med = swap_smallest_from_a_to_b(RML_rem, ids_bas, ids_med) # from medial to basal

    ids_bas, ids_med = smooth_edge(RML_rem, ids_bas, ids_med)
    ids_lat, ids_med = smooth_edge(RML_rem, ids_lat, ids_med)

    ids_bas, ids_med = clip_RML_base_to_med(RML_rem, ids_bas, ids_med, cent_RL_fiss)
    ids_bas, ids_med = swap_smallest_from_a_to_b(RML_rem, ids_bas, ids_med) # from basal to medial

    RML_bas = RML_rem.extract_cells(ids_bas).extract_largest()
    RML_lat = RML_rem.extract_cells(ids_lat).extract_largest()
    RML_med = RML_rem.extract_cells(ids_med).extract_largest()

    RML_lat_edge = edge_points(RML_lat)
    RML_med_edge = edge_points(RML_med)
    RML_base_edge = edge_points(RML_bas)

    edge_RM_base_lat = annotate_edge(RML_lat_edge, RML_base_edge, average_length/2)
    edge_RM_base_med = annotate_edge(RML_med_edge, RML_base_edge, average_length/2)
    edge_RM_anterior = annotate_edge(RML_lat_edge, RML_med_edge, average_length/2)
    edge_RM_H_fiss_lat = annotate_edge(RU_M_fis_edge, RML_lat_edge, average_length/2)
    edge_RM_H_fiss_med = annotate_edge(RU_M_fis_edge, RML_med_edge, average_length/2)

    debug = False
    if debug:
        p = pv.Plotter()
        pv.global_theme.window_size = [1000, 1000]
        p.add_mesh(RML_rem, color="black", style='wireframe')
        p.add_mesh(RUL_med, color="orange")
        p.add_mesh(RUL_lat, color="magenta")
        p.add_mesh(RLL_rem, color="blue", style='wireframe')
        p.add_mesh(RLL_med, color="red")
        p.add_mesh(RLL_lat, color="green")
        p.add_mesh(RLL_bas, color="orange")
        p.add_mesh(RML_med, color="green")
        p.add_mesh(RML_lat, color="red")
        p.add_mesh(RML_bas, color="gray")
        #p.add_points(ant_base_med_RLL, render_points_as_spheres=True, color='green', point_size = 20)
        #p.add_points(base_pt, render_points_as_spheres=True, color='pink', point_size = 20)
        p.camera_position = "xy"
        p.show()
        p.close()
    # -------------------------------------------------------------------------- end RML

    summed_RUL = RUL_med.n_cells + RUL_lat.n_cells
    summed_RML = RML_med.n_cells + RML_lat.n_cells + RML_bas.n_cells
    summed_RLL = RLL_med.n_cells + RLL_lat.n_cells + RLL_bas.n_cells
    err_RUL = 100*(1 - summed_RUL/RUL_rem.n_cells)
    err_RML = 100*(1 - summed_RML/RML_rem.n_cells)
    err_RLL = 100*(1 - summed_RLL/RLL_rem.n_cells)
    print (" RUL err: %8.2f %%  %6d  %6d" % (err_RUL, summed_RUL, RUL_rem.n_cells))
    print (" RML err: %8.2f %%  %6d  %6d" % (err_RML, summed_RML, RML_rem.n_cells))
    print (" RLL err: %8.2f %%  %6d  %6d" % (err_RLL, summed_RLL, RLL_rem.n_cells))
    print (" Number of edge points: ")
    print (" RUL --  ant_post (%d), Hfiss_lat (%d), Hfiss_med (%d) "
           % ( len(edge_RU_ant_post), len(edge_RU_M_fiss_lat), len(edge_RU_M_fiss_med)))
    print (" RML --  base_lat (%d), base_med (%d), anterior (%d) "
           % ( len(edge_RM_base_lat), len(edge_RM_base_med), len(edge_RM_anterior)))
    print (" RLL --  base_lat (%d), base_med (%d), fiss_base (%d), posterior (%d),"
           % (len(edge_RL_base_lat), len(edge_RL_base_med), len(edge_RL_fiss_base), len(edge_RL_posterior)))
    print ("         Hfiss_lat (%d), Hfiss_med (%d), Ofiss_lat (%d), Ofiss_med (%d)"
           % (len(edge_RL_U_fiss_lat), len(edge_RL_U_fiss_med), len(edge_RL_M_fiss_lat), len(edge_RL_M_fiss_med) ))
    
    # ---- OUTPUT THE ANNOTATED GROUPS ----

    # surfaces
    RUL_med_points = get_connected_points(RUL_med)
    RUL_lat_points = get_connected_points(RUL_lat)
    RML_med_points = get_connected_points(RML_med)
    RML_lat_points = get_connected_points(RML_lat)
    RLL_med_points = get_connected_points(RLL_med)
    RLL_lat_points = get_connected_points(RLL_lat)
    RLL_bas_points = get_connected_points(RLL_bas)
    output_mesh(os.path.join(output_directory, 'RUL_medial'), RUL_med_points)
    output_mesh(os.path.join(output_directory, 'RUL_lateral'), RUL_lat_points)
    output_mesh(os.path.join(output_directory, 'RML_medial'), RML_med_points)
    output_mesh(os.path.join(output_directory, 'RML_lateral'), RML_lat_points)
    output_mesh(os.path.join(output_directory, 'RML_base'), RML_bas.points[:])
    output_mesh(os.path.join(output_directory, 'RLL_medial'), RLL_med_points)
    output_mesh(os.path.join(output_directory, 'RLL_lateral'), RLL_lat_points)
    output_mesh(os.path.join(output_directory, 'RLL_base'), RLL_bas_points)

    pv.save_meshio(os.path.join(output_directory, 'RUL_medial.stl'), RUL_med, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RUL_medial.stl'), RUL_med, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RUL_lateral.stl'), RUL_lat, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RML_medial.stl'), RML_med, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RML_lateral.stl'), RML_lat, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RML_base.stl'), RML_bas, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RLL_medial.stl'), RLL_med, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RLL_lateral.stl'), RLL_lat, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RLL_base.stl'), RLL_bas, binary=True)

    # fissure surfaces
    output_mesh(os.path.join(output_directory, 'R_oblique_fissure'), R_oblique_fiss.points[:])
    output_mesh(os.path.join(output_directory, 'R_horiz_fissure'), R_horiz_fiss.points[:])
    pv.save_meshio(os.path.join(output_directory, 'RU_M_fissure.stl'), RU_M_fiss_surface, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RU_L_fissure.stl'), RU_L_fiss_surface, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'RL_M_fissure.stl'), RL_M_fiss_surface, binary=True)

    # edges around the RLL
    output_mesh(os.path.join(output_directory, 'edge_RL_base_lat'), edge_RL_base_lat)
    output_mesh(os.path.join(output_directory, 'edge_RL_base_med'), edge_RL_base_med)
    output_mesh(os.path.join(output_directory, 'edge_RL_fiss_base'), edge_RL_fiss_base)
    output_mesh(os.path.join(output_directory, 'edge_RL_posterior'), edge_RL_posterior)
    output_mesh(os.path.join(output_directory, 'edge_RL_H_fiss_lat'), edge_RL_U_fiss_lat)
    output_mesh(os.path.join(output_directory, 'edge_RL_H_fiss_med'), edge_RL_U_fiss_med)
    output_mesh(os.path.join(output_directory, 'edge_RM_O_fiss_lat'), edge_RL_M_fiss_lat)
    output_mesh(os.path.join(output_directory, 'edge_RM_O_fiss_med'), edge_RL_M_fiss_med)

    # edges around the RML
    output_mesh(os.path.join(output_directory, 'edge_RM_base_lat'), edge_RM_base_lat)
    output_mesh(os.path.join(output_directory, 'edge_RM_base_med'), edge_RM_base_med)
    output_mesh(os.path.join(output_directory, 'edge_RM_anterior'), edge_RM_anterior)

    # edges for the RUL
    output_mesh(os.path.join(output_directory, 'edge_RU_ant_post'), edge_RU_ant_post)
    output_mesh(os.path.join(output_directory, 'edge_RU_M_fiss_lat'), edge_RU_M_fiss_lat)
    output_mesh(os.path.join(output_directory, 'edge_RU_M_fiss_med'), edge_RU_M_fiss_med)

    debug = False
    if debug:
        p = pv.Plotter()
        pv.global_theme.window_size = [1600, 1600]
        p.add_mesh(RUL_med, color="orange")
        p.add_mesh(RUL_lat, color="magenta")
        p.add_mesh(RLL_med, color="blue")
        p.add_mesh(RLL_lat, color="green")
        p.add_mesh(RLL_bas, color="gold")
        p.add_mesh(RML_med, color="yellow")
        p.add_mesh(RML_lat, color="red")
        p.add_mesh(RML_bas, color="gray")
        p.add_points(apex_RUL, render_points_as_spheres=True, color='green', point_size = 20)
        p.add_points(anterior_point, render_points_as_spheres=True, color='orange', point_size = 20)
        p.add_points(inflexion_point, render_points_as_spheres=True, color='yellow', point_size = 20)
        p.show()
    
#############################################################################################################

def left_surface_data(output_directory, phenotype, dist_lim, sharp):

    LLL_surf = pv.read(os.path.join(output_directory, 'LLL_surf.stl'))
    LUL_surf = pv.read(os.path.join(output_directory, 'LUL_surf.stl'))

    #---- separate the fissures from the lobes

    _ = LLL_surf.compute_implicit_distance(LUL_surf, inplace=True)
    _ = LUL_surf.compute_implicit_distance(LLL_surf, inplace=True)
    LLL_surf['implicit_distance'] = np.abs(LLL_surf['implicit_distance'])
    LUL_surf['implicit_distance'] = np.abs(LUL_surf['implicit_distance'])

    LLL_fiss_surface = LLL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)
    LLL_rem = LLL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False)

    LUL_fiss_surface = LUL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=True)
    LUL_rem = LUL_surf.clip_scalar(scalars='implicit_distance', value = dist_lim, invert=False)
    LUL_rem = LUL_rem.smooth_taubin(n_iter=10, edge_angle=60, boundary_smoothing=False,
                                    feature_smoothing=True, non_manifold_smoothing=True)
    LLL_rem = LLL_rem.smooth_taubin(n_iter=10, edge_angle=60, boundary_smoothing=False,
                                    feature_smoothing=True, non_manifold_smoothing=True)
    
    #---- separate the LLL base from lateral and medial
    #     use ray tracing to find the cells that a vector from the cell to a point
    #     beneath the base will pass through only once (i.e. through the cell itself)
    
    bounds_min = np.min(LLL_rem.points, axis=0)
    bounds_max = np.max(LLL_rem.points, axis=0)
    bounds_centre = 0.5 * (bounds_max + bounds_min)
    base_centre = np.array([bounds_centre[0], bounds_centre[1], bounds_min[2] - 10.0])

    LLL_rem.point_data.clear()
    cell_centers = LLL_rem.cell_centers()   # geometric centre of each cell
    ids_base = np.empty((0), dtype=int)     # cell ids for 'base' and (below) 'other'
    ids_other = np.empty((0), dtype=int)
    
    LLL_rem.compute_normals(inplace=True, cell_normals=True, point_normals=False)
    cell_normals = LLL_rem['Normals'] # should be the normals at vertices
    cell_centers = LLL_rem.cell_centers()

    x_min = bounds_min[0] + (bounds_max[0]-bounds_min[0])*0.25
    x_max = bounds_max[0]
    xdrn = np.array([1.0, 0.0, 0.0])
    x_angle_limit = math.pi * 90.0/180.0
    ids_bas = []
    ids_lat = []
    ids_med = []

    for cell in range(LLL_rem.n_cells):
        #  using ray tracing to identify cells on base
        point, index = LLL_rem.ray_trace(cell_centers.points[cell], base_centre)
        if len(index) == 0 or (len(index) == 1 and index[0] == cell):
            ids_bas.append(cell)
        else:
            x = cell_centers.points[cell,0]
            x_angle = math.acos(np.dot(xdrn, cell_normals[cell]))
            if x < x_max and x > x_min and abs(x_angle) < x_angle_limit:
                ids_lat.append(cell)
            else:
                ids_med.append(cell)

    # swap isolated cells between the two groups
    ids_bas, ids_med = smooth_edge(LLL_rem, ids_bas, ids_med) # only swap into the base (pick up new ones)
    ids_lat, ids_bas = smooth_edge(LLL_rem, ids_lat, ids_bas) # only swap out to the lateral
    #ids_bas, ids_med = check_for_swaps(LLL_rem, ids_bas, ids_med) # only swap into the base
    #ids_bas, ids_med = check_for_swaps(LLL_rem, ids_bas, ids_med) # only swap into the base (pick up new ones)
    #ids_lat, ids_bas = check_for_swaps(LLL_rem, ids_lat, ids_bas) # only swap out to the lateral

    # create new surface entities using all cells from the id lists
    LLL_bas = LLL_rem.extract_cells(ids_bas).extract_largest().extract_surface()
    
    #---- identify landmark locations on LLL
    base_centre = LLL_bas.center_of_mass()
    av_base_normal = LLL_bas.compute_normals(cell_normals=True)['Normals'].mean(axis=0)
    av_base_normal /= np.linalg.norm(av_base_normal)
    ind_max = np.argmax(LLL_bas.points, axis=0)
    ind_min = np.argmin(LLL_bas.points, axis=0)
    post_base_LLL = np.array(LLL_bas.points[ind_max[1]]) # most posterior point
    ind_max = np.argmax(LLL_fiss_surface.points, axis=0)
    post_apex_LLL = np.array(LLL_fiss_surface.points[ind_max[2]])

    # get the cells that fall on each side of a dividing plane, based
    # on location of apical, and posterior LLL apical points
    COM_LLL = LLL_rem.center_of_mass()
    normal, d = plane_from_3pts(post_apex_LLL, post_base_LLL, COM_LLL)
    y_lim = 0.8*bounds_min[1] + 0.2*bounds_max[1] #COM_LLL[1]
    z_lim = post_base_LLL[2] - 30
    ids_med, ids_lat = swap_cells_using_plane(LLL_rem, cell_centers, ids_med, ids_lat,
                                              normal, d, y_lim, z_lim, True, True)
    ids_lat, ids_med = smooth_edge(LLL_rem, ids_lat, ids_med) # to lateral
    ids_lat, ids_med = swap_smallest_from_a_to_b(LLL_rem, ids_lat, ids_med) # from lateral to medial
    ids_bas, ids_med = swap_smallest_from_a_to_b(LLL_rem, ids_bas, ids_med) # from lateral to medial

    # create new surface entities using all cells from the id lists
    LLL_med = LLL_rem.extract_cells(ids_med).extract_largest()
    LLL_lat = LLL_rem.extract_cells(ids_lat).extract_largest()

    ind_max = np.argmax(LUL_rem.points, axis=0)
    ind_min = np.argmin(LUL_rem.points, axis=0)
    apex_LUL = np.array(LUL_rem.points[ind_max[2]])

    LLL_bas_edge = edge_points(LLL_bas)
    LLL_lat_edge = edge_points(LLL_lat)
    LLL_med_edge = edge_points(LLL_med)
    LLL_fis_edge = edge_points(LLL_fiss_surface)

    mean_surface_area = LLL_rem.area/LLL_rem.n_cells
    average_length = math.sqrt((4*mean_surface_area)/math.sqrt(3))
    
    edge_LL_base_lat = annotate_edge(LLL_bas_edge, LLL_lat_edge, average_length/2)
    edge_LL_base_med = annotate_edge(LLL_bas_edge, LLL_med_edge, average_length/2)
    edge_LL_fiss_base = annotate_edge(LLL_bas_edge, LLL_fis_edge, average_length/2)
    edge_LL_posterior = annotate_edge(LLL_lat_edge, LLL_med_edge, average_length/2)
    edge_LL_fiss_lat = annotate_edge(LLL_fis_edge, LLL_lat_edge, average_length/2)
    edge_LL_fiss_med = annotate_edge(LLL_fis_edge, LLL_med_edge, average_length/2)

    # check for phenotype
    if len(edge_LL_fiss_base) == 0: 
        phenotype_1 = 1
    else:
        phenotype_1 = 0

    if phenotype_1 != 1:
        point_LL_fiss_lat = get_closest_between(edge_LL_base_lat, edge_LL_fiss_base)
        point_LL_fiss_med = get_closest_between(edge_LL_base_med, edge_LL_fiss_base)

    #---- separate the LLL base from lateral and medial
    #     use normal directions to estimate the allocation of cells
    
    LUL_rem.compute_normals(inplace=True, cell_normals=False, point_normals=True)
    LUL_normals = LUL_rem['Normals'] # should be the normals at vertices
    bounds_min = np.min(LUL_rem.points, axis=0)
    bounds_max = np.max(LUL_rem.points, axis=0)
    bounds_centre = 0.5 * (bounds_max + bounds_min)
    bounds_range = (bounds_max - bounds_min)
    
    LUL_rem.point_data.clear()
    LUL_rem.compute_normals(inplace=True, cell_normals=True, point_normals=False)
    cell_normals = LUL_rem['Normals'] # should be the normals at vertices
    cell_centers = LUL_rem.cell_centers()

    # find the 'knee' on anterior, taken as the most apical point in the anterior region
    # (defined as y < y_lim) that will be classified as medial (using ray tracing)
    COM_LUL = LUL_rem.center_of_mass()
    knee = np.array([-1.0e8,-1.0e8,-1.0e8])
    centre_LUL = np.array([bounds_min[0], COM_LUL[1], COM_LUL[2]])
    for cell in range(LUL_rem.n_cells):
        x = cell_centers.points[cell,0]
        y = cell_centers.points[cell,1]
        z = cell_centers.points[cell,2]
        if x < COM_LUL[0] and y < y_lim:
            #  using ray tracing to identify cells on LUL medial
            point, index = LUL_rem.ray_trace(cell_centers.points[cell], centre_LUL)
            if len(index) == 0 or (len(index) == 1 and index[0] == cell):
                if z > knee[2]:
                    knee = cell_centers.points[cell]


    lat_centre1 = np.array([bounds_max[0] + 0.5*bounds_range[0], bounds_min[1]
                           - 0.3*bounds_range[1], bounds_centre[2] + 0.5*bounds_range[2]])
    lat_centre2 = np.array([bounds_max[0] + 0.5*bounds_range[0], bounds_min[1]
                           - 0.6*bounds_range[1], bounds_centre[2]])

    xyz_range = bounds_max - bounds_min
    xyz_factor = [0.25, 0.35, 0.2]
    x_min = bounds_min[0] + xyz_range[0] * xyz_factor[0]
    x_max = bounds_max[0]
    xdrn = np.array([1.0, 0.0, 0.0])
    y_min = bounds_min[1] # i.e. no limit
    y_max = bounds_min[1] + xyz_range[1] * xyz_factor[1]
    ydrn = np.array([0.0, -1.0, 0.0])
    z_min = bounds_min[2] # i.e. no limit
    if phenotype_1 != 1:
        z_max = max(point_LL_fiss_lat[2], point_LL_fiss_med[2]) + 10.0
    zdrn = np.array([0.0, 0.0, -1.0])
    x_angle_limit = math.pi * 90.0/180.0
    y_angle_limit = math.pi * 60.0/180.0
    z_angle_limit = math.pi * 60.0/180.0
    # angle 60 for z (base) used to not pick up too many false positives
    angle_limit = [value * math.pi/180.0 for value in [60.0, 60.0, 60.0]]
    ant_apical = np.array((1.0e8, 1.0e8,-1.0e8))
    c_angle_min = -1.0e8
    ids_bas, ids_lat, ids_med = [], [], []
    ids_med2 = []

    for cell in range(LUL_rem.n_cells):
        x = cell_centers.points[cell,0]
        y = cell_centers.points[cell,1]
        z = cell_centers.points[cell,2]
        x_angle = math.acos(np.dot(xdrn, cell_normals[cell]))
        y_angle = math.acos(np.dot(ydrn, cell_normals[cell]))
        z_angle = math.acos(np.dot(zdrn, cell_normals[cell]))
        #z_angle = math.acos(np.dot(av_base_normal, cell_normals[cell]))
        bdrn = (base_centre - cell_centers.points[cell]) # use the angle to the base_centre
        bdrn = bdrn/np.linalg.norm(bdrn)
        b_angle = math.acos(np.dot(bdrn, cell_normals[cell]))
        cdrn = cell_centers.points[cell] - bounds_centre
        cdrn[2] = 0.0
        cdrn = cdrn/np.linalg.norm(cdrn)
        c_angle = math.acos(np.dot(cdrn, xdrn))
        # note that ray tracing not used for LUL base because it is only a small region, so
        # not much of the region 'above' it is shadowed (too many false positives)
        
        if x < x_max and x > x_min and abs(x_angle) < angle_limit[0]:
            ids_lat.append(cell)
        else:
            if z > bounds_centre[2] + 0.2*bounds_range[2]:
                point, index = LUL_rem.ray_trace(cell_centers.points[cell], lat_centre1)
                if len(index) == 0 or (len(index) == 1 and index[0] == cell):
                    ids_lat.append(cell)
                else:
                    ids_med.append(cell)
            else:
                if y < y_max and 0.0 < y_angle < angle_limit[1]:
                    ids_lat.append(cell)
                    if z > ant_apical[2]:
                        ant_apical[2] = z
                    if c_angle > c_angle_min:
                        c_angle_min = c_angle
                        ant_apical[0], ant_apical[1] = x, y
                elif phenotype_1 != 1: # only when LUL_base exists
                    dist_lfm = math.dist(point_LL_fiss_med, cell_centers.points[cell])
                    # accept as basal if not too far from the med-fiss point or lower than its
                    # z height, plus meets angle criteria
                    if (dist_lfm < 10.0 or z < z_max) and (abs(b_angle) < math.radians(60.0)
                                                           or abs(z_angle) < angle_limit[2]):
                        ids_bas.append(cell)
                    else:
                        ids_med.append(cell)
                else:
                    ids_med.append(cell)

    ant_apical_cell = get_closest_point_to_group(ant_apical, cell_centers.points)
    ant_apical = cell_centers.points[ant_apical_cell]
    
    # get the cells that fall on each side of a dividing plane, based
    # on location of apical and posterior LLL apical points
    COM_LUL = LUL_rem.center_of_mass()
    normal, d = plane_from_3pts(apex_LUL, post_apex_LLL, COM_LUL)
    y_lim = apex_LUL[1]
    z_lim = post_apex_LLL[2] - 30
    ids_med, ids_lat = swap_cells_using_plane(LUL_rem, cell_centers, ids_med, ids_lat,
                                              normal, d, y_lim, z_lim, True, True)

    normal, d = plane_from_3pts(apex_LUL, ant_apical, COM_LUL)
    y_lim = apex_LUL[1]
    ids_med, ids_lat = swap_cells_using_plane(LUL_rem, cell_centers, ids_med, ids_lat,
                                              normal, d, y_lim, z_lim, False, True)

    ids_lat, ids_med = swap_smallest_between(LUL_rem, ids_lat, ids_med) # between lateral and medial
    ids_bas, ids_med = swap_smallest_between(LUL_rem, ids_bas, ids_med) # between basal and lateral

    # clip apical 'knotch' using a clipping plane
    frac_ant_pnt = 0.7  # location of the anterior point along z-axis of LUL
    frac_inf_pnt = 0.9  # location of the inflexion point along z-axis of LUL
    anterior_point, inflexion_point = get_anterior_inflexion(COM_LUL, LUL_rem, ids_lat, frac_ant_pnt, frac_inf_pnt)

    # clip between inflexion point and anterior
    normal, d = plane_from_3pts(COM_LUL, inflexion_point, anterior_point)
    y_max = apex_LUL[1]
    y_min = anterior_point[1]
    z_min = anterior_point[2]
    ids_med, ids_lat = swap_cells_from_a_to_b(cell_centers, ids_med, ids_lat, normal, d, y_min, y_max, z_min, False)

    # clip between inflexion point and apical point
    normal, d = plane_from_3pts(COM_LUL, apex_LUL, inflexion_point)
    y_min = inflexion_point[1]
    z_min = inflexion_point[2]
    ids_lat, ids_med = swap_cells_from_a_to_b(cell_centers, ids_lat, ids_med, normal, d, y_min, y_max, z_min, True)

    # smooth the edge between adjacent surfaces
    if phenotype_1 != 1:
        ids_bas, ids_med = smooth_edge(LUL_rem, ids_bas, ids_med)
        ids_bas, ids_lat = smooth_edge(LUL_rem, ids_bas, ids_lat)
    ids_lat, ids_med = smooth_edge(LUL_rem, ids_lat, ids_med)

    # create new surface entities using all cells from the id lists
    if phenotype_1 != 1:
        LUL_bas = LUL_rem.extract_cells(ids_bas).extract_largest()
    LUL_lat = LUL_rem.extract_cells(ids_lat).extract_largest()
    LUL_med = LUL_rem.extract_cells(ids_med).extract_largest()

    if phenotype_1 != 1:
        LUL_bas_edge = edge_points(LUL_bas)
    LUL_lat_edge = edge_points(LUL_lat)
    LUL_med_edge = edge_points(LUL_med)
    LUL_fis_edge = edge_points(LUL_fiss_surface)
    
    if phenotype_1 != 1:
        edge_LU_base_lat = annotate_edge(LUL_bas_edge, LUL_lat_edge, average_length/2)
        edge_LU_base_med = annotate_edge(LUL_bas_edge, LUL_med_edge, average_length/2)
        edge_LU_fiss_base = annotate_edge(LUL_bas_edge, LUL_fis_edge, average_length/2)
    edge_LU_post_ant = annotate_edge(LUL_lat_edge, LUL_med_edge, average_length/2)

    
    # -------------------------------------------------------------------------- end LUL

    if phenotype_1 != 1:
        summed_LUL = LUL_med.n_cells + LUL_lat.n_cells + LUL_bas.n_cells
    else:
        summed_LUL = LUL_med.n_cells + LUL_lat.n_cells
    summed_LLL = LLL_med.n_cells + LLL_lat.n_cells + LLL_bas.n_cells
    err_LUL = 100*(1 - summed_LUL/LUL_rem.n_cells)
    err_LLL = 100*(1 - summed_LLL/LLL_rem.n_cells)
    print (" LUL err: %8.2f %%  %6d  %6d" % (err_LUL, summed_LUL, LUL_rem.n_cells))
    print (" LLL err: %8.2f %%  %6d  %6d" % (err_LLL, summed_LLL, LLL_rem.n_cells))
    print (" Number of edge points: ")
    if phenotype_1 != 1:
        print (" LUL --  ant_post (%d), base_lat (%d), base_med (%d) "
               % ( len(edge_LU_post_ant), len(edge_LU_base_lat), len(edge_LU_base_med)))
    else:
        print (" LUL --  ant_post (%d) " % ( len(edge_LU_post_ant)))
    print (" LLL --  base_lat (%d), base_med (%d), fiss_base (%d), posterior (%d),"
           % (len(edge_LL_base_lat), len(edge_LL_base_med), len(edge_LL_fiss_base), len(edge_LL_posterior)))
    print ("         fiss_lat (%d), fiss_med (%d)"
           % (len(edge_LL_fiss_lat), len(edge_LL_fiss_med) ))
    
    # ---- OUTPUT THE ANNOTATED GROUPS ----

    # surfaces
    LUL_med_points = get_connected_points(LUL_med)
    LUL_lat_points = get_connected_points(LUL_lat)
    if phenotype_1 != 1: LUL_bas_points = get_connected_points(LUL_bas)
    LLL_med_points = get_connected_points(LLL_med)
    LLL_lat_points = get_connected_points(LLL_lat)
    LLL_bas_points = get_connected_points(LLL_bas)
    output_mesh(os.path.join(output_directory, 'LUL_medial'), LUL_med_points)
    output_mesh(os.path.join(output_directory, 'LUL_lateral'), LUL_lat_points)
    if phenotype_1 != 1: output_mesh(os.path.join(output_directory, 'LUL_base'), LUL_bas_points)
    output_mesh(os.path.join(output_directory, 'LLL_medial'), LLL_med_points)
    output_mesh(os.path.join(output_directory, 'LLL_lateral'), LLL_lat_points)
    output_mesh(os.path.join(output_directory, 'LLL_base'), LLL_bas_points)

    pv.save_meshio(os.path.join(output_directory, 'LUL_medial.stl'), LUL_med, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'LUL_lateral.stl'), LUL_lat, binary=True)
    if phenotype_1 != 1: pv.save_meshio(os.path.join(output_directory, 'LUL_base.stl'), LUL_bas, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'LLL_medial.stl'), LLL_med, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'LLL_lateral.stl'), LLL_lat, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'LLL_base.stl'), LLL_bas, binary=True)

    # fissure surfaces
    output_mesh(os.path.join(output_directory, 'LUL_fissure'), LUL_fiss_surface.points[:])
    output_mesh(os.path.join(output_directory, 'LLL_fissure'), LLL_fiss_surface.points[:])
    pv.save_meshio(os.path.join(output_directory, 'LUL_fissure.stl'), LUL_fiss_surface, binary=True)
    pv.save_meshio(os.path.join(output_directory, 'LLL_fissure.stl'), LLL_fiss_surface, binary=True)

    # edges around the LUL
    if phenotype_1 != 1: output_mesh(os.path.join(output_directory, 'edge_LU_base_lat'), edge_LU_base_lat)
    if phenotype_1 != 1: output_mesh(os.path.join(output_directory, 'edge_LU_base_med'), edge_LU_base_med)
    if phenotype_1 != 1: output_mesh(os.path.join(output_directory, 'edge_LU_fiss_base'), edge_LU_fiss_base)
    output_mesh(os.path.join(output_directory, 'edge_LU_post_ant'), edge_LU_post_ant)

    # edges around the LLL
    output_mesh(os.path.join(output_directory, 'edge_LL_base_lat'), edge_LL_base_lat)
    output_mesh(os.path.join(output_directory, 'edge_LL_base_med'), edge_LL_base_med)
    output_mesh(os.path.join(output_directory, 'edge_LL_fiss_base'), edge_LL_fiss_base)
    output_mesh(os.path.join(output_directory, 'edge_LL_posterior'), edge_LL_posterior)
    output_mesh(os.path.join(output_directory, 'edge_LL_fiss_lat'), edge_LL_fiss_lat)
    output_mesh(os.path.join(output_directory, 'edge_LL_fiss_med'), edge_LL_fiss_med)

    debug = False
    if debug:
        p = pv.Plotter()
        pv.global_theme.window_size = [1200, 1200]
        p.add_mesh(LLL_rem, color="orange", style='wireframe')
        if phenotype_1 != 1: # read groups for the LUL base
            p.add_mesh(LUL_bas, color="orange")
            #p.add_mesh(LUL_neck, color="white")
            #p.add_mesh(LUL_2, color="blue")
            #p.add_mesh(LUL_3, color="green", style='wireframe')
            #if len(edge_LU_base_lat) > 0: p.add_points(edge_LU_base_lat, render_points_as_spheres=True, color='white', point_size=10)
            #if len(edge_LU_base_med) > 0: p.add_points(edge_LU_base_med, render_points_as_spheres=True, color='green', point_size=10)
            #if len(edge_LU_fiss_base) > 0: p.add_points(edge_LU_fiss_base, render_points_as_spheres=True, color='blue', point_size=10)
        p.add_mesh(LUL_lat, color="red")
        p.add_mesh(LUL_med, color="magenta")
        p.add_mesh(LLL_bas, color="cyan")
        p.add_mesh(LLL_lat, color="yellow")
        p.add_mesh(LLL_med, color="green")
        #p.add_points(point_LL_fiss_lat, render_points_as_spheres=True, color='green', point_size = 20)
        #p.add_points(point_LL_fiss_med, render_points_as_spheres=True, color='red', point_size = 20)
        #p.add_points(anterior_point, render_points_as_spheres=True, color='orange', point_size = 20)
        #p.add_points(inflexion_point, render_points_as_spheres=True, color='yellow', point_size = 20)
        #p.add_points(apex_LUL, render_points_as_spheres=True, color='green', point_size = 20)
        p.show()
    
#############################################################################################################

# identify edge points that are closer to a sequentially 'distant' edge point than to neighbours. From this,
# get the list of cell ids (for cells located closer to the medial side) that use these points. Remove these
# cells from ids_bas and add to ids_med. Doesn't matter that the removed/added cells are not contiguous: later
# operations that swap all smallest connected groups from the base to medial will take care of this. 

def clip_RML_base_to_med(RML_rem, ids_bas, ids_med, cent_RL_fiss):

    list_lines = get_list_of_lines_ids(RML_rem, ids_bas) # all lines that are on the RML base

    # find the edge lines
    edge = {}
    num_edges = 0
    edge_points = {}
    keep_lines = []

    for cell_id in range(RML_rem.n_cells):
        if cell_id in ids_bas: # we only want edges that are in the base
            pt_ids = RML_rem.cell_point_ids(cell_id)
            pt_ids.sort()
            lines = [f"{pt_ids[0]}_{pt_ids[1]}", f"{pt_ids[0]}_{pt_ids[2]}", f"{pt_ids[1]}_{pt_ids[2]}"]
            for i, line in enumerate(lines):
                count_lines = list_lines.count(line)
                if count_lines == 1 and line not in keep_lines:
                    keep_lines.append(line)
                    start, end = line.split('_')
                    edge[num_edges] = [int(start), int(end), cell_id, None, None]
                    num_edges += 1
    # identify the neighbouring edges/nodes
    for key in edge:
        start, end = edge[key][:2]  # Get start and end values of current edge
        edge[key][3] = next((k for k in edge if k != key and (start == edge[k][0] or start == edge[k][1])), None)
        edge[key][4] = next((k for k in edge if k != key and (end == edge[k][0] or end == edge[k][1])), None)
        if start == edge[edge[key][3]][0]:
            edge_points[start] = [end, edge[edge[key][3]][1], edge[key][2]]
        elif start == edge[edge[key][3]][1]:
            edge_points[start] = [end, edge[edge[key][3]][0], edge[key][2]]
        if end == edge[edge[key][4]][0]:
            edge_points[end] = [start, edge[edge[key][4]][1], edge[key][2]]
        elif end == edge[edge[key][4]][1]:
            edge_points[end] = [start, edge[edge[key][4]][0], edge[key][2]]

    num_keys = len(edge_points)
    ordered_edge = {}
    ordpnts = []
    ordcell = []
    point = min(edge_points.keys())   # start from the edge point of smallest value
    point0, point1 = edge_points[point][0], edge_points[point][1]
    ordpnts.extend([point0, point, point1])  
    ordcell.extend([edge_points[point0][2], edge_points[point][2], edge_points[point1][2]]) # the edge cell ids
    
    while len(ordpnts) < num_keys:
        for k in (edge_points[point1][0], edge_points[point1][1]):
            if k != point:
                point2 = k
                ordpnts.append(point2)  
                ordcell.append(edge_points[point2][2])
                break  # Exit the loop after finding the value
        point = point1
        point1 = point2

    # get a list of point coordinates that are 'too close'
    return_coords = np.empty((0,3), dtype=float)
    bridge_elements = np.empty((0,3), dtype=int) # stores lines joining the points
    list_cells = []
    point_ids = []

    for i, point in enumerate(ordpnts):
        point_coords = RML_rem.points[point]
        if point_coords[0] > cent_RL_fiss[0]:
            close_list = []
            for k in range(5):
                p = ordpnts[i-k]
                close_list.append(p)
                if i+k < len(ordpnts):
                    p = ordpnts[i+k]
                else:
                    p = ordpnts[k]
                close_list.append(p)

            coord_list = np.empty((0,3), dtype=float)
            point_list = [] # a list of points to check the distance to
            for j, p in enumerate(ordpnts):
                if p not in close_list:
                    coord_list = np.vstack((coord_list, RML_rem.points[p]))
                    point_list.append(p)
            point_coord = RML_rem.points[point]
            distances = [np.linalg.norm(np.array(point_coord) - np.array(p)) for p in coord_list]
            shortest = min(distances)
            edge_length = np.linalg.norm(np.array(point_coord) - np.array(RML_rem.points[close_list[3]]))
            if 0.1 < shortest < edge_length:
                j = np.argmin(distances) # the index in 'distances' array of the closest point
                point_ids.append(point)
                point_ids.append(point_list[j])

    for cell_id in ids_bas:
        for point in point_ids:
            if point in RML_rem.cell_point_ids(cell_id):
                list_cells.append(cell_id)
            
    list_cells = list(set(list_cells)) # remove duplicates
    for cell in list_cells:
        ids_med.append(cell)
        ids_bas.remove(cell)

    return ids_bas, ids_med 
                
#############################################################################################################

def py_mesh_and_repair(mesh, resolution):

    # function to improve mesh quality while also reducing number of faces
    mesh, target_len = fix_mesh(mesh, detail=resolution)

    # access vertices and faces from mesh
    pyv_vertices = mesh.vertices
    mesh_faces = mesh.faces
    pyv_faces = np.zeros((mesh.num_faces*4), dtype=int)
    ncount = 0
    for i in range(mesh.num_faces):
        pyv_faces[ncount] = 3
        ncount += 1
        for j in range(3):
            pyv_faces[ncount] = mesh_faces[i,j]
            ncount += 1

    # make a polydata surface for visualising with PyVista
    surf = pv.PolyData(pyv_vertices, pyv_faces)
    
    return surf, pyv_vertices, target_len

    
#############################################################################################################

def process_data_from_stl(input_data, output_directory, lung, lobes, phenotype,
                          img_scale, dist_fiss, annotate, read_if_exist, sharp):
    
    resolution = "low"

    LUL_stl = os.path.join(input_data, lobes['LUL'] + '.stl')
    LLL_stl = os.path.join(input_data, lobes['LLL'] + '.stl')
    RUL_stl = os.path.join(input_data, lobes['RUL'] + '.stl')
    RML_stl = os.path.join(input_data, lobes['RML'] + '.stl')
    RLL_stl = os.path.join(input_data, lobes['RLL'] + '.stl')
    
    make_files = True
    # read stl files for surfaces, extract surface and edge data, write to files
    if lung == 'Left':
        path_to_file = os.path.join(output_directory, 'LLL_surf.stl')
        if read_if_exist == True and os.path.exists(path_to_file) == True:
            make_files = False
        if make_files:
            # remake the files, or they don't exist
            print(' Make lobe surface files')
            LLL_mesh = pymesh.load_mesh(LLL_stl)
            LLL_surf, LLL_vert, target_LLL = py_mesh_and_repair(LLL_mesh, resolution) # returns unscaled data
            
            LUL_mesh = pymesh.load_mesh(LUL_stl)
            LUL_surf, LUL_vert, target_LUL = py_mesh_and_repair(LUL_mesh, resolution) # returns unscaled data
            
            opv = open(os.path.join(output_directory, 'left_stl_volumes.txt'), 'w+')
            opv.write(" STL volume LLL = %8.2f mL " % (LLL_surf.volume*1.0e-3 *
                                                       img_scale[0]*img_scale[1]*img_scale[2]))
            opv.write(" STL volume LUL = %8.2f mL " % (LUL_surf.volume*1.0e-3 *
                                                       img_scale[0]*img_scale[1]*img_scale[2]))
            opv.close()
    
            # write out the surfaces so that this does not have to be repeated
            
            LLL_surf = LLL_surf.fill_holes(hole_size=100)
            LUL_surf = LUL_surf.fill_holes(hole_size=100)

            LLL_surf = LLL_surf.scale(img_scale, inplace=True)
            LUL_surf = LUL_surf.scale(img_scale, inplace=True)

            # smoothing causes issue with identifying the fissures. sharp fingers
            LLL_surf.smooth(inplace=True,n_iter=50, feature_smoothing=True, edge_angle=60, feature_angle=60)
            LUL_surf.smooth(inplace=True,n_iter=50, feature_smoothing=True, edge_angle=60, feature_angle=60)

            LLL_surf.save(os.path.join(output_directory, 'LLL_surf.ply'), binary=False)
            LUL_surf.save(os.path.join(output_directory, 'LUL_surf.ply'), binary=False)
    
            pv.save_meshio(os.path.join(output_directory, 'LLL_surf.stl'), LLL_surf, binary=True)
            pv.save_meshio(os.path.join(output_directory, 'LUL_surf.stl'), LUL_surf, binary=True)
            # end 'if make_files'
            debug = False
            if debug:
                p = pv.Plotter()
                pv.global_theme.window_size = [1200, 1200]
                p.add_mesh(LLL_surf, color="black", style='wireframe')
                p.add_mesh(LLL_surf, color="orange")
                p.add_mesh(LUL_surf, color="black", style='wireframe')
                p.add_mesh(LUL_surf, color="green")
                p.show()

        if annotate:
            print ('Annotating surfaces')
            # process the lobe surfaces to identify fissures, outer surfaces
            #   this identifies and writes out the R_oblique_fiss, R_horiz_fiss, RM_oblique_fiss, RM_horiz_fiss
            #   and RLL, RML, RUL remaining after subtracting the fissures
            left_surface_data(output_directory, phenotype, dist_fiss, sharp) 

    elif lung == 'Right':
        path_to_file = os.path.join(output_directory, 'RLL_surf.stl')
        if read_if_exist == True and os.path.exists(path_to_file) == True:
            make_files = False
        if make_files:
            # remake the files, or they don't exist
            print(' Make lobe surface files')
            RLL_mesh = pymesh.load_mesh(RLL_stl)
            RLL_surf, RLL_vert, target_RLL = py_mesh_and_repair(RLL_mesh, resolution) # returns unscaled data

            RML_mesh = pymesh.load_mesh(RML_stl)
            RML_surf, RML_vert, target_RML = py_mesh_and_repair(RML_mesh, resolution) # returns unscaled data

            RUL_mesh = pymesh.load_mesh(RUL_stl)
            RUL_surf, RUL_vert, target_RUL = py_mesh_and_repair(RUL_mesh, resolution) # returns unscaled data

            opv = open(os.path.join(output_directory, 'right_stl_volumes.txt'), 'w+')
            opv.write(" STL volume RLL = %8.2f mL " % (RLL_surf.volume*1.0e-3 *
                                                       img_scale[0]*img_scale[1]*img_scale[2]))
            opv.write(" STL volume RML = %8.2f mL " % (RML_surf.volume*1.0e-3 *
                                                       img_scale[0]*img_scale[1]*img_scale[2]))
            opv.write(" STL volume RUL = %8.2f mL " % (RUL_surf.volume*1.0e-3 *
                                                       img_scale[0]*img_scale[1]*img_scale[2]))
            opv.close()

            # write out the surfaces so that this does not have to be repeated
            
            RLL_surf = RLL_surf.fill_holes(hole_size=100)
            RML_surf = RML_surf.fill_holes(hole_size=100)
            RUL_surf = RUL_surf.fill_holes(hole_size=100)

            RLL_surf = RLL_surf.scale(img_scale, inplace=True)
            RML_surf = RML_surf.scale(img_scale, inplace=True)
            RUL_surf = RUL_surf.scale(img_scale, inplace=True)

            # smoothing causes issue with identifying the fissures. sharp fingers
            RLL_surf.smooth(inplace=True,n_iter=50, feature_smoothing=True, edge_angle=60, feature_angle=60)
            RML_surf.smooth(inplace=True,n_iter=50, feature_smoothing=True, edge_angle=60, feature_angle=60)
            RUL_surf.smooth(inplace=True,n_iter=50, feature_smoothing=True, edge_angle=60, feature_angle=60)

            RLL_surf.save(os.path.join(output_directory, 'RLL_surf.ply'), binary=False)
            RML_surf.save(os.path.join(output_directory, 'RML_surf.ply'), binary=False)
            RUL_surf.save(os.path.join(output_directory, 'RUL_surf.ply'), binary=False)
    
            pv.save_meshio(os.path.join(output_directory, 'RLL_surf.stl'), RLL_surf, binary=True)
            pv.save_meshio(os.path.join(output_directory, 'RML_surf.stl'), RML_surf, binary=True)
            pv.save_meshio(os.path.join(output_directory, 'RUL_surf.stl'), RUL_surf, binary=True)
            # end 'if make_files'

        if annotate:
            # process the lobe surfaces to identify fissures, outer surfaces
            #   this identifies and writes out the R_oblique_fiss, R_horiz_fiss, RM_oblique_fiss, RM_horiz_fiss
            #   and RLL, RML, RUL remaining after subtracting the fissures
            right_surface_data(output_directory, phenotype, dist_fiss, sharp) 






