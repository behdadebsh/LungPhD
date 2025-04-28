import trimesh
from trimesh.ray.ray_triangle import ray_triangle_id
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import json
import os
from natsort import natsorted
from collections import defaultdict
import copy
import matplotlib.cm as cm
import random
import pyvista as pv
import matplotlib
matplotlib.use('TkAgg')


def load_json_curves(json_files, directory):
    # Extract control points and positions into a dictionary
    node_dict = {}
    current_id = 1
    # Loop through all JSON files in the directory
    for json_count, filename in enumerate(json_files):
        file_path = os.path.join(directory, filename)

        # Load JSON file
        with open(file_path, "r") as file:
            data = json.load(file)
            temp_dict = {}
            temp_count = 0
            for markup in data.get("markups", []):
                for point in markup.get("controlPoints", []):
                    temp_count += 1
                    temp_dict[temp_count] = {
                        "branch": json_count,
                        "coords": [point["position"][0], point["position"][1], point["position"][2]]
                    }

            # Downsample, ensuring first and last points are kept
            resampled_points = downsample_points(temp_dict, step=8)
            renumbered_points, current_id = renumber_points(resampled_points, current_id)
            node_dict.update(renumbered_points)
        json_count += 1

    return node_dict


def downsample_points(points, step):
    """Downsample while keeping the first and last point."""
    sorted_keys = sorted(points.keys(), key=int)  # Ensure numeric order
    if len(sorted_keys) <= step:
        downsampled = {sorted_keys[0]: points[sorted_keys[0]]}
        downsampled[sorted_keys[-1]] = points[sorted_keys[-1]]  # Keep last point
    else:
        downsampled = {sorted_keys[0]: points[sorted_keys[0]]}  # Keep first point
        downsampled.update({k: points[k] for k in sorted_keys[step:-1:step]})  # Downsample middle
        downsampled[sorted_keys[-1]] = points[sorted_keys[-1]]  # Keep last point

    return downsampled


def renumber_points(points, start_id):
    """Renumber points sequentially from start_id and return new start_id."""
    renumbered = {str(i + start_id): coords for i, (_, coords) in enumerate(points.items())}
    new_start_id = start_id + len(renumbered)  # Update counter for next file
    return renumbered, new_start_id


def identify_branch_points(node_dict):
    """Identify branch points by finding duplicate nodes"""
    # Finding duplicate coordinates
    coordinate_map = defaultdict(list)
    branch_map = defaultdict(list)
    first_occurrence = {}
    for node, data in node_dict.items():
        coord = tuple(data["coords"])
        branch = data["branch"]

        if coord not in first_occurrence:
            first_occurrence[coord] = node  # Mark first occurrence

        coordinate_map[coord].append((node, branch))
        branch_map[branch].append((node, coord))
    # Add duplicate info to the first occurrence node
    for coord, node_list in coordinate_map.items():
        if len(node_list) > 1:  # Only process duplicates
            first_node = first_occurrence[coord]
            duplicates = [{"node": n, "branch": c} for n, c in node_list if n != first_node]
            node_dict[first_node]["duplicates"] = duplicates  # Add duplicate info

            # Delete all other duplicate nodes
            for node, _ in node_list:
                if node != first_node:
                    del node_dict[node]

    # Re-number the nodes after deleting the duplicates
    node_dict_renumber = {}
    for i, occurrence in enumerate(first_occurrence.items()):
        node_dict_renumber[str(i + 1)] = node_dict[occurrence[1]]

    return node_dict_renumber


def define_parent_child(node_dict):
    node_keys = list(node_dict.keys())  # Get a list of node names
    for i, node_name in enumerate(node_keys):
        current_branch = node_dict[node_name]["branch"]

        if i == 0:  # First node is inlet - no parent
            node_dict[node_name]["parent"] = []
            node_dict[node_name]["child"] = [node_keys[i + 1]]
        elif i == len(node_keys) - 1:  # Last node is always terminal - no child
            node_dict[node_name]["parent"] = node_keys[i - 1]
            node_dict[node_name]["child"] = []
        else:  # Other nodes
            # If next node is part of current branch
            if node_dict[node_keys[i + 1]]["branch"] == current_branch:
                node_dict[node_name]["child"] = [node_keys[i + 1]]
            else:
                if "child" not in node_dict[node_name]:
                    if "duplicates" in node_dict[node_name]:  # Is a branch point
                        duplicates = node_dict[node_name]["duplicates"]
                        branch_values = [node["branch"] for node in duplicates]

                        first_node = []
                        for branch in branch_values:
                            nodes_in_branch = {node_id for node_id, data in node_dict.items() if data["branch"] == branch}
                            temp_node = min(nodes_in_branch, key=int)
                            node_dict[temp_node]["parent"] = node_name  # Set parent node for branching child node (as cannot be tracked later)
                            first_node.append(temp_node)

                        node_dict[node_name]["child"] = first_node  # Set child nodes for this branch point

                    else:  # Terminal node
                        node_dict[node_name]["child"] = []

            # If parent node not already defined
            if "parent" not in node_dict[node_name]:
                node_dict[node_name]["parent"] = node_keys[i - 1]

    return node_dict


def find_closest_pair(numbers):
    min_distance = float('inf')
    indices = None

    # Compare every pair of numbers
    for i in range(len(numbers)):
        for j in range(i + 1, len(numbers)):
            dist = abs(numbers[i] - numbers[j])  # Use absolute difference for numbers
            if dist < min_distance:
                min_distance = dist
                indices = [i, j]

    remaining_indices = [i for i in range(len(numbers)) if i not in indices]
    return indices, remaining_indices


def farthest_index(nums):
    # Calculate the absolute differences between each number and the others
    distances = [
        (abs(nums[0] - nums[1]), abs(nums[0] - nums[2])),  # Distances for nums[0]
        (abs(nums[1] - nums[0]), abs(nums[1] - nums[2])),  # Distances for nums[1]
        (abs(nums[2] - nums[0]), abs(nums[2] - nums[1])),  # Distances for nums[2]
    ]
    # Calculate the total distance for each number
    total_distances = [sum(dist) for dist in distances]

    return total_distances.index(max(total_distances))


def fix_multiple_branches(node_dict):
    node_keys = list(node_dict.keys())  # Get a list of node names
    for node_name in node_keys:
        if len(node_dict[node_name]["child"]) == 3:
            # Fix for trifurcation - move one branch up to parent node
            # Need to add check if parent has bifurcation already
            duplicates = node_dict[node_name]["duplicates"]
            branches = [node["branch"] for node in duplicates]
            parent = node_dict[node_name]["parent"]
            child = node_dict[node_name]["child"]

            max_ind = farthest_index(branches)
            child_move = child[max_ind]

            # Update parent and child for affected nodes
            child_rest = child[:max_ind] + child[max_ind + 1:]
            if len(node_dict[parent]["child"]) > 1:
                stop = 0
                i = 0
                while stop == 0:
                    if len(node_dict[child_rest[i]]["child"]) == 1:
                        new_parent = child_rest[i]
                        node_dict[child_move]["parent"] = new_parent
                        node_dict[new_parent]["child"].append(child_move)
                        stop = 1
                    i += 1
                    if i == len(child_rest):
                        stop = 1
            else:
                node_dict[child_move]["parent"] = parent
                node_dict[parent]["child"].append(child_move)

            node_dict[node_name]["child"] = child_rest

        elif len(node_dict[node_name]["child"]) == 4:
            duplicates = node_dict[node_name]["duplicates"]
            branches = [node["branch"] for node in duplicates]
            parent = node_dict[node_name]["parent"]
            child = node_dict[node_name]["child"]

            ind_keep, ind_move = find_closest_pair(branches)
            child_keep = [child[i] for i in ind_keep]
            # child_move = [child[i] for i in ind_move]

            # Update parent and child for affected nodes
            for temp_ind in ind_move:
                child_rest = child[:temp_ind] + child[temp_ind + 1:]
                child_move = child[temp_ind]
                if len(node_dict[parent]["child"]) > 1:
                    stop = 0
                    i = 0
                    while stop == 0:
                        if len(node_dict[child_rest[i]]["child"]) == 1:
                            new_parent = child_rest[i]
                            node_dict[child_move]["parent"] = new_parent
                            node_dict[new_parent]["child"].append(child_move)
                            stop = 1
                        i += 1
                        if i == len(child_rest):
                            stop = 1
                else:
                    node_dict[child_move]["parent"] = parent
                    node_dict[parent]["child"].append(child_move)

            node_dict[node_name]["child"] = child_keep

    return node_dict


def writeExNodeFile(node_dict, node_file):
    """
    Write out an ex Node file with the given coords.
    :param node_file: Filename to write to.
    :param nodes: List of coordinate lists.
    :return: None
    """
    with open(node_file, 'w') as f:
        f.write(' Group name : Upper_tree\n')
        f.write(' #Fields=1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 0\n')
        f.write('   y.  Value index= 2, #Derivatives= 0\n')
        f.write('   z.  Value index= 3, #Derivatives= 0\n')

        for node, keys in node_dict.items():
            f.write(' Node:         %s\n' % node)
            f.write('    %3.5f\n' % float(keys["coords"][0]))
            f.write('    %3.5f\n' % float(keys["coords"][1]))
            f.write('    %3.5f\n' % float(keys["coords"][2]))


def writeExElemFile(node_dict, elem_file):
    """
        Write out an ex Node file with the given coords.
        :param filename: Filename to write to.
        :param coords: List of coordinate lists.
        :return: None
    """
    with open(elem_file, 'w') as f:
        f.write(' Group name: Upper_tree\n')
        f.write(' Shape.  Dimension=1\n')
        f.write(' #Scale factor sets=1\n')
        f.write('  l.lagrange, #Scale factors= 2\n')
        f.write(' #Nodes= 2\n')
        f.write(' #Fields=1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  l.Lagrange, no modify, standard node based.\n')
        f.write('     #Nodes= 2\n')
        f.write('      1.  #Values=1\n')
        f.write('       Value indices:   1\n')
        f.write('       Scale factor indices:   1\n')
        f.write('      2.  #Values=1\n')
        f.write('       Value indices:   1\n')
        f.write('       Scale factor indices:   2\n')
        f.write('   y.  l.Lagrange, no modify, standard node based.\n')
        f.write('     #Nodes= 2\n')
        f.write('     1.  #Values=1\n')
        f.write('       Value indices:   1\n')
        f.write('       Scale factor indices:   1\n')
        f.write('     2.  #Values=1\n')
        f.write('       Value indices:   1\n')
        f.write('       Scale factor indices:   2\n')
        f.write('   z.  l.Lagrange, no modify, standard node based.\n')
        f.write('     #Nodes= 2\n')
        f.write('     1.  #Values=1\n')
        f.write('       Value indices:   1\n')
        f.write('       Scale factor indices:   1\n')
        f.write('     2.  #Values=1\n')
        f.write('       Value indices:   1\n')
        f.write('       Scale factor indices:   2\n')

        count = 0
        for node, keys in node_dict.items():
            child_all = node_dict[node]["child"]
            parent_node = node
            for child_node in child_all:
                count += 1
                f.write(' Element:         %.d 0 0\n' % count)
                f.write('   Nodes:\n')
                f.write('       %s         %s\n' % (parent_node, child_node))
                f.write('   Scale factors:\n')
                f.write('     0.10000E+01   0.10000E+01\n')


def get_meta_from_json(file_path):
    """
    Reads a Slicer markups JSON file and returns the first 'orientation' matrix found.
    """
    with open(file_path, 'r') as file:
        data = json.load(file)

    for markup in data.get("markups", []):
        for cp in markup.get("controlPoints", []):
            if "orientation" in cp:
                orientation = cp["orientation"]
                return np.array([orientation[0:3], orientation[3:6], orientation[6:9]])

    return None


def label_from_vector(vec):
    axis_names = [
        ("Left", "Right"),
        ("Posterior", "Anterior"),
        ("Superior", "Inferior")
    ]
    idx = np.argmax(np.abs(vec))
    sign = np.sign(vec[idx])
    return axis_names[idx][0] if sign > 0 else axis_names[idx][1], axis_names[idx][1] if sign > 0 else axis_names[idx][0], f"{axis_names[idx][1]}-{axis_names[idx][0]}"


def visualize_node_graph(data, orientation_matrix=None, arrow_scale=0.1, color_by_branch=True,
                         ply_path=None, terminals=None):
    """
    Visualizes the vascular network in 3D and overlays anatomical orientation and terminal nodes with PLY mesh.

    Args:
        data: The vascular tree dictionary
        orientation_matrix: 3x3 orientation matrix (columns = i, j, k directions)
        arrow_scale: Arrow scale for flow direction
        color_by_branch: Whether to color by branch ID
        ply_path: Optional path to a closed surface mesh (.ply)
        terminals: Optional list of terminal nodes [(node_id, element)] to visualize
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("3D Vascular Tree with Orientation")

    # Optional: branch color map
    branch_colors = {}
    cmap = cm.get_cmap('tab20', 100)

    def get_branch_color(branch_id):
        if branch_id not in branch_colors:
            branch_colors[branch_id] = cmap(random.randint(0, 99))
        return branch_colors[branch_id]

    # Plot vascular structure
    for node_id, node_data in data.items():
        start = np.array(node_data['coords'])
        for child_id in node_data.get('child', []):
            if child_id not in data:
                continue
            child = np.array(data[child_id]['coords'])
            branch_id = node_data['branch'] if color_by_branch else 0
            color = get_branch_color(branch_id)

            ax.plot([start[0], child[0]], [start[1], child[1]], [start[2], child[2]], color=color, linewidth=1.5)

            # Directional arrow
            direction = child - start
            norm = np.linalg.norm(direction)
            if norm > 0:
                direction /= norm
                ax.quiver(start[0], start[1], start[2],
                          direction[0], direction[1], direction[2],
                          length=arrow_scale * norm, color=color, arrow_length_ratio=0.3, linewidth=0.5)

    # Plot anatomical axes if provided
    if orientation_matrix is not None:
        center = np.mean(np.array([v['coords'] for v in data.values()]), axis=0)
        orientation = np.array(orientation_matrix)
        scale = 50
        for i in range(3):
            vec = orientation[:, i]
            label, opp_label, label_range = label_from_vector(vec)

            ax.quiver(center[0], center[1], center[2],
                      -vec[0], -vec[1], vec[2],
                      length=scale, color='red', linewidth=2, arrow_length_ratio=0.1)
            ax.text(center[0] + vec[0] * scale,
                    center[1] + vec[1] * scale,
                    center[2] + vec[2] * scale,
                    f"{label}\n({label_range})", color='black', fontsize=10)
            ax.text(center[0] - vec[0] * scale,
                    center[1] - vec[1] * scale,
                    center[2] - vec[2] * scale,
                    f"{opp_label}", color='gray', fontsize=9)

    # Overlay terminal node inside/outside status
    if ply_path and terminals:
        mesh = pv.read(ply_path)
        points = np.array([data[node_id]['coords'] for node_id, _ in terminals])
        point_cloud = pv.PolyData(points)
        enclosed = mesh.select_enclosed_points(point_cloud, tolerance=1e-6, check_surface=True)
        inside_flags = enclosed['SelectedPoints'][:len(points)].astype(bool)

        inside_pts = points[inside_flags]
        outside_pts = points[~inside_flags]

        # Show in plot
        if len(inside_pts) > 0:
            ax.scatter(inside_pts[:, 0], inside_pts[:, 1], inside_pts[:, 2],
                       color='blue', label='Terminals Inside', s=20, alpha=1)
        if len(outside_pts) > 0:
            ax.scatter(outside_pts[:, 0], outside_pts[:, 1], outside_pts[:, 2],
                       color='red', label='Terminals Outside', s=20, alpha=1)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.tight_layout()
    plt.show()


def visualize_mesh_with_tree_and_terminals(tree_dict, terminals, ply_paths):
    """
    Visualize multiple PLY meshes with vascular tree and terminal points overlaid,
    including labels for terminal node IDs, and assign distinct colors to each mesh.
    """
    # Create point and line sets for the vascular tree
    points = []
    lines = []
    id_map = {}  # Map node_id to index in points list
    idx = 0

    for node_id, node in tree_dict.items():
        coord = node['coords']
        points.append(coord)
        id_map[node_id] = idx
        idx += 1

    for node_id, node in tree_dict.items():
        for child_id in node.get('child', []):
            if child_id in tree_dict:
                lines.append([2, id_map[node_id], id_map[child_id]])  # 2-point line

    # Build PyVista Line Mesh
    points_np = np.array(points)
    lines_np = np.hstack(lines)
    tree = pv.PolyData()
    tree.points = points_np
    tree.lines = lines_np

    # Terminal points and labels
    term_coords = np.array([tree_dict[node_id]['coords'] for node_id, _ in terminals])
    term_labels = [str(node_id) for node_id, _ in terminals]

    # Setup plotter
    plotter = pv.Plotter()

    # Define a list of distinct colors
    color_list = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan', 'orange', 'purple', 'brown', 'pink']

    # Add each mesh from the list with a distinct color
    for i, path in enumerate(ply_paths):
        mesh = pv.read(path)
        color = color_list[i % len(color_list)]  # Cycle through colors if more meshes than colors
        plotter.add_mesh(mesh, color=color, opacity=0.3, show_edges=False, label=f"Mesh {i+1}")

    plotter.add_mesh(tree, color="black", line_width=1.5, label="Vascular Tree")
    plotter.add_points(term_coords, color="red", point_size=10, render_points_as_spheres=True, label="Terminals")

    # Add labels to terminal nodes with enhanced visibility settings
    plotter.add_point_labels(
        term_coords,
        term_labels,
        font_size=12,
        text_color="black",
        point_color="red",
        point_size=10,
        render_points_as_spheres=True,
        always_visible=True,
        pickable=False
    )

    plotter.add_legend()
    plotter.show()


def find_terminals(tree_dict):
    """
        Finds terminal nodes in a vascular tree dictionary.
        Returns a list of tuples: (terminal_node_id, parent_element)
        where parent_element is (parent_node_id, terminal_node_id)
        """
    terminals = []

    for node_id, node_data in tree_dict.items():
        # No children = terminal node
        if 'child' not in node_data or not node_data['child']:
            parent = node_data.get('parent')
            if parent:
                terminals.append((node_id, (parent, node_id)))
            else:
                terminals.append((node_id, None))  # root with no children?
    return terminals


def find_inlet_nodes(tree_dict):
    """
    Finds inlet nodes in a vascular tree dictionary.
    Returns a list of tuples: (inlet_node_id, first_element)
    where first_element is (inlet_node_id, first_child_id)
    """
    inlets = []

    for node_id, node_data in tree_dict.items():
        parent = node_data.get('parent')
        if not parent:  # no parent or empty list = inlet
            child = node_data.get('child')
            if child:
                inlets.append((node_id, (node_id, child[0])))
            else:
                inlets.append((node_id, None))  # orphan node with no parent
    if len(inlets) > 1:
        print("WARNING! More that one inlet detected. Things might crash or not work onwards.")
    return inlets


def check_terminals_inside_mesh(terminals, tree_dict, ply_path):
    """
    Checks which terminal nodes are inside a mesh using ray casting.
    """
    # Load the mesh
    mesh = trimesh.load(ply_path, process=True)
    triangles = mesh.triangles  # (N, 3, 3)

    # Extract terminal coordinates
    points = np.array([tree_dict[node_id]['coords'] for node_id, _ in terminals])

    # Define ray directions (e.g., positive X direction)
    ray_directions = np.tile([1.0, 0.0, 0.0], (len(points), 1))

    # Perform ray-triangle intersection
    index_tri, index_ray, _ = ray_triangle_id(
        triangles=triangles,
        ray_origins=points,
        ray_directions=ray_directions,
        multiple_hits=True
    )

    # Count intersections per ray
    hit_counts = np.bincount(index_ray, minlength=len(points))

    # Determine if points are inside based on odd number of intersections
    inside_flags = hit_counts % 2 == 1

    # Compile results
    results = []
    for i, (node_id, _) in enumerate(terminals):
        results.append((node_id, points[i], bool(inside_flags[i])))

    return results


def check_unassigned_terminals(terminals, terminal_to_lobe):
    """
    Identifies terminal nodes not assigned to any lung lobe.

    Args:
        terminals (list): List of tuples, each containing a node ID and its connection, e.g., [('23', ('22', '23')), ...].
        terminal_to_lobe (dict): Dictionary mapping terminal node IDs to lobe names.

    Returns:
        list: Terminal node IDs not assigned to any lobe.
    """
    # Extract node IDs from the terminals list
    terminal_ids = {node_id for node_id, _ in terminals}
    # Get the set of assigned terminal IDs
    assigned_ids = set(terminal_to_lobe.keys())
    # Determine unassigned terminal IDs
    unassigned = terminal_ids - assigned_ids
    return list(unassigned)


def build_parent_map(tree_dict):
    """
    Build a mapping from each node to its parent.
    """
    parent_map = {}
    for parent_id, node in tree_dict.items():
        for child_id in node.get('child', []):
            parent_map[child_id] = parent_id
    return parent_map


def get_ancestor_path(node_id, parent_map):
    """
    Retrieve the path from the given node up to the root.
    """
    path = [node_id]
    while node_id in parent_map:
        node_id = parent_map[node_id]
        path.append(node_id)
    return path


def find_lcp_for_lobes(tree_dict, terminal_to_lobe):
    """
    For each lung lobe, find the lowest common parent of its terminal nodes.
    """
    parent_map = build_parent_map(tree_dict)

    # Organize terminal nodes by lobe
    lobe_terminals = {}
    for terminal_id, lobe in terminal_to_lobe.items():
        lobe_terminals.setdefault(lobe, []).append(terminal_id)

    lobe_lca = {}
    for lobe, terminals in lobe_terminals.items():
        # Get ancestor paths for all terminals in the lobe
        ancestor_paths = [get_ancestor_path(tid, parent_map) for tid in terminals]

        # Find common ancestors by intersecting ancestor sets
        common_ancestors = set(ancestor_paths[0])
        for path in ancestor_paths[1:]:
            common_ancestors &= set(path)

        if not common_ancestors:
            lobe_lca[lobe] = None
            continue

        # Determine the deepest common ancestor (i.e., the one closest to the terminals)
        # This is the common ancestor with the minimal depth in any of the paths
        # Since paths are from node to root, the first common node encountered is the LCA
        for node in ancestor_paths[0]:
            if node in common_ancestors:
                lobe_lca[lobe] = node
                break
    return lobe_lca


def find_feeding_elements(tree_dict, lobe_lcp_mapping):
    """
    For each lung lobe, find the parent-child connection (edge) that feeds all its terminal nodes.

    Args:
        tree_dict (dict): Dictionary representing the vascular tree. Each key is a node ID, and its value is a dictionary with keys like 'coords' and 'child'.
        lobe_lca_mapping (dict): Dictionary mapping each lobe name to its lowest common ancestor (LCA) node ID.

    Returns:
        dict: Mapping from lobe names to tuples representing the feeding edge (parent_id, lca_id).
    """
    # Build a parent map for quick lookup
    parent_map = {}
    for parent_id, node in tree_dict.items():
        for child_id in node.get('child', []):
            parent_map[child_id] = parent_id

    # Determine the feeding edge for each lobe
    lobe_feeding_edges = {}
    for lobe, lca_id in lobe_lcp_mapping.items():
        parent_id = parent_map.get(lca_id)
        if parent_id:
            lobe_feeding_edges[lobe] = (parent_id, lca_id)
        else:
            # If the LCA is the root node (has no parent), the feeding edge is undefined
            lobe_feeding_edges[lobe] = None

    return lobe_feeding_edges


def assign_rul_subsegments(centroids, orientation_matrix):
    # Project centroids into anatomical space
    projected = centroids @ orientation_matrix
    z_vals = projected[:, 2]  # Superior-Inferior
    y_vals = projected[:, 1]  # Anterior-Posterior

    indices = np.arange(len(centroids))
    assignments = {}

    # Step 1: Identify Apical clusters
    apical_indices = indices[np.argsort(-z_vals)[:2]]
    if y_vals[apical_indices[0]] < y_vals[apical_indices[1]]:
        assignments[apical_indices[0]] = "Apical Superior"
        assignments[apical_indices[1]] = "Apical Inferior"
    else:
        assignments[apical_indices[1]] = "Apical Superior"
        assignments[apical_indices[0]] = "Apical Inferior"

    # Remaining indices
    remaining_indices = list(set(indices) - set(apical_indices))

    # Step 2: Identify Anterior clusters
    anterior_candidates = sorted(remaining_indices, key=lambda i: y_vals[i])[:2]
    if z_vals[anterior_candidates[0]] > z_vals[anterior_candidates[1]]:
        assignments[anterior_candidates[0]] = "Anterior Superior"
        assignments[anterior_candidates[1]] = "Anterior Inferior"
    else:
        assignments[anterior_candidates[1]] = "Anterior Superior"
        assignments[anterior_candidates[0]] = "Anterior Inferior"

    # Remaining indices
    remaining_indices = list(set(remaining_indices) - set(anterior_candidates))

    # Step 3: Identify Posterior clusters
    if z_vals[remaining_indices[0]] > z_vals[remaining_indices[1]]:
        assignments[remaining_indices[0]] = "Posterior Superior"
        assignments[remaining_indices[1]] = "Posterior Inferior"
    else:
        assignments[remaining_indices[1]] = "Posterior Superior"
        assignments[remaining_indices[0]] = "Posterior Inferior"

    return assignments



root_folder = '/hpc/bsha219/lung/Data'
study = 'CTEPH'
protocol = 'Pre'
subject_name = 'Alfred12'
subject_path = os.path.join(root_folder, study, subject_name, protocol)
ply_directory = os.path.join(root_folder, study, subject_name, protocol, 'Lobe', 'Lobe_mesh')
RUL_mesh = pv.read(os.path.join(ply_directory, 'RUL_surf.ply'))
RML_mesh = pv.read(os.path.join(ply_directory, 'RML_surf.ply'))
RLL_mesh = pv.read(os.path.join(ply_directory, 'RLL_surf.ply'))
LUL_mesh = pv.read(os.path.join(ply_directory, 'LUL_surf.ply'))
LLL_mesh = pv.read(os.path.join(ply_directory, 'LLL_surf.ply'))
lobe_mesh_paths = {
    "RUL": os.path.join(ply_directory, 'RUL_surf.ply'),
    "RML": os.path.join(ply_directory, 'RML_surf.ply'),
    "RLL": os.path.join(ply_directory, 'RLL_surf.ply'),
    "LUL": os.path.join(ply_directory, 'LUL_surf.ply'),
    "LLL": os.path.join(ply_directory, 'LLL_surf.ply')
}

# Directory containing JSON files
json_directory = os.path.join(
    root_folder, study, subject_name, protocol, 'Vessel',
    'Alfred12_Pre_VMTK')  # Change this to your directory path
vmtk_files = natsorted([f for f in os.listdir(json_directory) if f.startswith("Centerline curve")])


centreline_dict = load_json_curves(vmtk_files, json_directory)

centreline_dict_renumber = copy.deepcopy(centreline_dict)
centreline_dict_renumber = identify_branch_points(centreline_dict_renumber)

centreline_dict_hierachy = copy.deepcopy(centreline_dict_renumber)
centreline_dict_hierachy = define_parent_child(centreline_dict_hierachy)

centreline_dict_final = copy.deepcopy(centreline_dict_hierachy)
centreline_dict_final = fix_multiple_branches(centreline_dict_final)
terminals = find_terminals(centreline_dict_final)

# Initialize a dictionary to hold the mapping
terminal_to_lobe = {}

# Iterate over each lobe and its corresponding mesh
for lobe_name, mesh_path in lobe_mesh_paths.items():
    # Check which terminals are inside the current lobe's mesh
    results = check_terminals_inside_mesh(terminals, centreline_dict_final, mesh_path)

    # For each terminal node, if it's inside the mesh, record the lobe
    for node_id, coords, is_inside in results:
        if is_inside:
            terminal_to_lobe[node_id] = lobe_name

unassigned = check_unassigned_terminals(terminals, terminal_to_lobe)  # find if any of the upper tree terminals is not assigned to a lobe
if len(unassigned) > 0: # this bit of the code tries to assign them to the nearest lobe (surface)
    unassigned_coords = np.array([centreline_dict_final[node]['coords'] for node in unassigned])
    for idx, coord in enumerate(unassigned_coords):
        min_distance = np.inf
        closest_lobe = None
        coord = np.array(coord).reshape(1, 3)  # Make it 2D array because find_closest_cell expects arrays

        for lobe_name, ply_path in lobe_mesh_paths.items():
            mesh = pv.read(ply_path)
            _, closest_point = mesh.find_closest_cell(coord, return_closest_point=True)

            # Calculate Euclidean distance
            distance = np.linalg.norm(coord - closest_point)

            if distance < min_distance:
                min_distance = distance
                closest_lobe = lobe_name

        terminal_to_lobe[unassigned[idx]] = closest_lobe

# lobe_lca_mapping = find_lcp_for_lobes(centreline_dict_final, terminal_to_lobe)   # option to find which lobe they belong to
# lobe_feeding_edges = find_feeding_elements(centreline_dict_final, lobe_lca_mapping)  # option to find joint parent elem feeding the lobe

to_export = False   # Option to export the CMISS tree
if to_export:
    output_exnode = os.path.join(subject_path, '/Vessel/VMTK/Alfred12_Pre_MPA.exnode')
    output_exelem = os.path.join(subject_path, '/Vessel/VMTK/Alfred12_Pre_MPA.exelem')
    writeExNodeFile(centreline_dict_final, output_exnode)
    writeExElemFile(centreline_dict_final, output_exelem)
# -----------------------------
# Load orientation and mesh data
# -----------------------------

to_visualise = True  # Option to visualise, Useful for debugging and verification
if to_visualise:
    directions = get_meta_from_json(os.path.join(subject_path, 'Vessel', 'Alfred12_Pre_VMTK', 'Centerline.json'))
    directions = np.array(directions)
    ply_paths = [lobe_mesh_paths['RUL'], lobe_mesh_paths['RML'], lobe_mesh_paths['RLL'],
                 lobe_mesh_paths['LUL'], lobe_mesh_paths['LLL']]
    visualize_mesh_with_tree_and_terminals(centreline_dict_final, terminals, ply_paths)
