import numpy as np
import pyvista as pv
import os


def calculate_centroid_from_exdata(file_path):
    """
    Reads a .exdata file and calculates the centroid of all 3D node coordinates.

    Args:
        file_path (str): Path to the .exdata file.

    Returns:
        np.ndarray: Centroid as a NumPy array [x, y, z].
    """
    with open(file_path, "r") as file:
        lines = file.readlines()

    coordinates = []
    i = 0
    while i < len(lines):
        if lines[i].strip().lower().startswith("node"):
            i += 1
            parts = lines[i].strip().split()
            if len(parts) == 3:
                try:
                    coord = list(map(float, parts))
                    coordinates.append(coord)
                except ValueError:
                    pass
        i += 1

    coords_array = np.array(coordinates)
    if coords_array.size == 0:
        raise ValueError(f"No coordinates found in file {file_path}")

    return coords_array.mean(axis=0)

def get_cluster_centroids(directory, cluster_range=range(1, 31)):
    """
    Calculates centroids for multiple .exdata cluster files.

    Args:
        directory (str): Directory where the .exdata files are stored.
        cluster_range (iterable): Range of cluster numbers to read.

    Returns:
        dict: Dictionary with cluster number as key and centroid as np.ndarray value.
    """
    centroids = {}
    for cluster_id in cluster_range:
        file_path = os.path.join(directory, f"Cl_{cluster_id}.exdata")
        try:
            centroid = calculate_centroid_from_exdata(file_path)
            centroids[cluster_id] = centroid
        except Exception as e:
            centroids[cluster_id] = f"Error: {e}"
    return centroids


def filter_centroids_by_lobe(centroid_dict, lobe):
    """
    Filters the centroid dictionary by the specified lobe label.

    Args:
        centroid_dict (dict): Dictionary of cluster_id -> centroid np.array
        lobe (str): One of 'RUL', 'RML', 'RLL', 'LLL', 'LUL'

    Returns:
        dict: Filtered dictionary with only the clusters belonging to that lobe.
    """
    lobe_clusters = {
        "RUL": list(range(1, 7)),
        "RML": list(range(7, 11)),
        "RLL": list(range(11, 17)),
        "LLL": list(range(17, 23)),
        "LUL": list(range(23, 31))
    }

    if lobe not in lobe_clusters:
        raise ValueError(f"Unknown lobe '{lobe}'. Choose from {list(lobe_clusters.keys())}")

    cluster_ids = lobe_clusters[lobe]
    return {cid: centroid_dict[cid] for cid in cluster_ids if cid in centroid_dict}


def annotate_RML(centroid_dict):
    rml_centroids = filter_centroids_by_lobe(centroid_dict, "RML")
    # Step 1: Find 2 posterior (smallest Y)
    sorted_by_y = sorted(rml_centroids.items(), key=lambda item: item[1][1])
    posterior = sorted_by_y[2:]
    anterior = sorted_by_y[:2]

    labels = {}

    # Step 2: Label anterior pair
    a_lateral, a_medial = sorted(anterior, key=lambda item: item[1][0])  # sort by X
    labels[a_lateral[0]] = "Lateral Anterior"
    labels[a_medial[0]] = "Medial Anterior"

    # Step 3: Label posterior pair
    p_lateral, p_medial = sorted(posterior, key=lambda item: item[1][0])  # sort by X
    labels[p_lateral[0]] = "Lateral Posterior"
    labels[p_medial[0]] = "Medial Posterior"

    return labels


def annotate_RUL(centroid_dict):
    rul_centroids = filter_centroids_by_lobe(centroid_dict, "RUL")

    # Step 1: Find 2 posterior (smallest Y)
    sorted_by_y = sorted(rul_centroids.items(), key=lambda item: item[1][1])
    posterior = sorted_by_y[4:]
    anterior = sorted_by_y[:2]
    apical = sorted_by_y[2:4]

    labels = {}

    # Step 2: Label anterior pair
    a_inferior, a_superior = sorted(anterior, key=lambda item: item[1][2])  # sort by Z
    labels[a_superior[0]] = "Anterior Superior"
    labels[a_inferior[0]] = "Anterior Inferior"

    # Step 3: Label posterior pair
    p_inferior, p_superior = sorted(posterior, key=lambda item: item[1][2])  # sort by Z
    labels[p_inferior[0]] = "Posterior Superior"
    labels[p_superior[0]] = "Posterior Inferior"

    # Step 3: Label apical pair
    ap_inferior, ap_superior = sorted(apical, key=lambda item: item[1][2])  # sort by Z
    labels[ap_inferior[0]] = "Apical Inferior"
    labels[ap_superior[0]] = "Apical Superior"

    return labels


def annotate_RLL(centroid_dict):
    rll_centroids = filter_centroids_by_lobe(centroid_dict, "RLL")
    labels = {}

    # Step 1: Posterior & Anterior (Y)
    sorted_by_y = sorted(rll_centroids.items(), key=lambda item: item[1][1])
    anterior = sorted_by_y[0]
    posterior = sorted_by_y[-1]
    labels[posterior[0]] = "Posterior"
    labels[anterior[0]] = "Anterior"

    # Step 2: Superior upper (Z)
    sorted_by_z = sorted(rll_centroids.items(), key=lambda item: item[1][2], reverse=True)
    superior_upper = sorted_by_z[0]
    labels[superior_upper[0]] = "Superior upper"

    # Step 3: Remaining clusters
    labeled_ids = {posterior[0], anterior[0], superior_upper[0]}
    remaining = {cid: pos for cid, pos in rll_centroids.items() if cid not in labeled_ids}

    # Step 4: Lateral = lowest X
    sorted_by_x = sorted(remaining.items(), key=lambda item: item[1][0])
    lateral = sorted_by_x[0]
    labels[lateral[0]] = "Lateral"

    # Step 5: Remaining two → decide Medial & Superior lower
    posterior_x = posterior[1][0]
    remaining_2 = sorted_by_x[1:]  # two left

    # Closest to posterior X → Medial
    dist_to_posterior = [abs(p[1][0] - posterior_x) for p in remaining_2]
    medial = remaining_2[np.argmin(dist_to_posterior)]
    superior_lower = remaining_2[np.argmax(dist_to_posterior)]

    labels[medial[0]] = "Medial"
    labels[superior_lower[0]] = "Superior lower"

    return labels


def annotate_LLL(centroid_dict):
    lll_centroids = filter_centroids_by_lobe(centroid_dict, "LLL")
    labels = {}

    # Step 1: Posterior (min Y)
    sorted_by_y = sorted(lll_centroids.items(), key=lambda item: item[1][1])
    posterior = sorted_by_y[-1]
    posterior_id, posterior_coords = posterior
    labels[posterior_id] = "Posterior"

    # Step 2: Superior upper (max Z)
    sorted_by_z = sorted(lll_centroids.items(), key=lambda item: item[1][2], reverse=True)
    superior_upper = sorted_by_z[0]
    labels[superior_upper[0]] = "Superior upper"

    # Step 3: Remaining 4
    labeled_ids = {posterior_id, superior_upper[0]}
    remaining = {cid: pos for cid, pos in lll_centroids.items() if cid not in labeled_ids}

    # Step 4: Lateral (cosine with +X axis, x > posterior x)
    posterior_x = posterior_coords[0]
    x_axis = np.array([1, 0, 0])
    candidate_vectors = {
        cid: pos - posterior_coords
        for cid, pos in remaining.items()
        if pos[0] > posterior_x
    }

    if not candidate_vectors:
        raise ValueError("No remaining centroid has X > posterior X for lateral direction.")

    cosine_scores = {
        cid: np.dot(vec, x_axis) / np.linalg.norm(vec)
        for cid, vec in candidate_vectors.items()
    }

    lateral_id = max(cosine_scores, key=cosine_scores.get)
    labels[lateral_id] = "Lateral"
    remaining.pop(lateral_id)

    # Step 5: Superior lower (max Z from remaining)
    sorted_remaining_by_z = sorted(remaining.items(), key=lambda item: item[1][2], reverse=True)
    superior_lower = sorted_remaining_by_z[0]
    labels[superior_lower[0]] = "Superior lower"
    remaining.pop(superior_lower[0])

    # Step 6: Last two: compare X
    remaining_items = list(remaining.items())
    if len(remaining_items) != 2:
        raise RuntimeError("Unexpected number of remaining centroids. Should be 2.")

    (id1, pos1), (id2, pos2) = remaining_items
    if pos1[0] > pos2[0]:
        labels[id1] = "Anterior"
        labels[id2] = "Medial"
    else:
        labels[id2] = "Anterior"
        labels[id1] = "Medial"

    return labels


def annotate_LUL(centroid_dict):
    lul_centroids = filter_centroids_by_lobe(centroid_dict, "LUL")
    labels = {}

    # Step 1: Sort all by Z
    sorted_by_z = sorted(lul_centroids.items(), key=lambda item: item[1][2])

    lingular_inferior = sorted_by_z[0]
    lingular_superior = sorted_by_z[1]
    apical_superior = sorted_by_z[-1]

    labels[lingular_inferior[0]] = "Lingular inferior"
    labels[lingular_superior[0]] = "Lingular superior"
    labels[apical_superior[0]] = "Apical superior"

    # Step 2: Remove these 3 from further processing
    excluded_ids = {lingular_inferior[0], lingular_superior[0], apical_superior[0]}
    remaining = {cid: pos for cid, pos in lul_centroids.items() if cid not in excluded_ids}

    # Step 3: Sort remaining by Y
    sorted_by_y = sorted(remaining.items(), key=lambda item: item[1][1])
    anterior1, anterior2 = sorted_by_y[:2]
    posterior1, posterior2 = sorted_by_y[-2:]
    maybe_apical_inferior = [item for item in sorted_by_y[2:-2]]

    # Step 4: Anterior assignments
    if anterior1[1][2] > anterior2[1][2]:
        labels[anterior1[0]] = "Anterior superior"
        labels[anterior2[0]] = "Anterior inferior"
    else:
        labels[anterior2[0]] = "Anterior superior"
        labels[anterior1[0]] = "Anterior inferior"

    # Step 5: Posterior assignments
    if posterior1[1][2] > posterior2[1][2]:
        labels[posterior1[0]] = "Posterior superior"
        labels[posterior2[0]] = "Posterior inferior"
    else:
        labels[posterior2[0]] = "Posterior superior"
        labels[posterior1[0]] = "Posterior inferior"

    # Step 6: Remaining cluster → Apical inferior
    if len(maybe_apical_inferior) != 1:
        raise RuntimeError("Expected exactly 1 cluster remaining for Apical inferior.")

    labels[maybe_apical_inferior[0][0]] = "Apical inferior"

    return labels


root_folder = '/hpc/bsha219/lung/Data'
study = 'CTEPH'
protocol = 'Pre'
subject_name = 'Alfred12'
subject_path = os.path.join(root_folder, study, subject_name, protocol)
ply_directory = os.path.join(root_folder, study, subject_name, protocol, 'Lobe', 'Lobe_mesh')

# Dictionary of PLY mesh paths
lobe_files = {
    "RUL": os.path.join(ply_directory, 'RUL_surf.ply'),
    "RML": os.path.join(ply_directory, 'RML_surf.ply'),
    "RLL": os.path.join(ply_directory, 'RLL_surf.ply'),
    "LUL": os.path.join(ply_directory, 'LUL_surf.ply'),
    "LLL": os.path.join(ply_directory, 'LLL_surf.ply')
}

# Cluster centroids (example)
centroid_dict = get_cluster_centroids(os.path.join(subject_path, 'Intensity_mapping'))

# Load meshes
lobe_meshes = {}
for lobe, path in lobe_files.items():
    if os.path.exists(path):
        lobe_meshes[lobe] = pv.read(path)

# Setup PyVista plotter
plotter = pv.Plotter()
colors = {
    "RUL": "lightblue",
    "RML": "lightgreen",
    "RLL": "lightcoral",
    "LUL": "plum",
    "LLL": "lightgray"
}

# Create plotter
plotter = pv.Plotter()

# Add all lobe meshes
for lobe, path in lobe_files.items():
    if os.path.exists(path):
        mesh = pv.read(path)
        plotter.add_mesh(mesh, color=colors.get(lobe, "white"), opacity=0.3, label=lobe)

# Prepare centroid coordinates and labels
points = np.array(list(centroid_dict.values()))
labels = [f"Cl {cid}" for cid in centroid_dict.keys()]

rml = annotate_RML(centroid_dict)
rul = annotate_RUL(centroid_dict)
rll = annotate_RLL(centroid_dict)
lll = annotate_LLL(centroid_dict)
lul = annotate_LUL(centroid_dict)


# Plot centroids and cluster numbers
plotter.add_points(points, color="red", point_size=10, render_points_as_spheres=True)
plotter.add_point_labels(points, labels, font_size=12,
        text_color="black",
        point_color="red",
        point_size=10,
        render_points_as_spheres=True,
        always_visible=True,
        pickable=False)

# Show legend and axes
plotter.add_legend()
plotter.show_axes()
plotter.show()
