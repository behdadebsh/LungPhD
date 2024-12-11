import numpy as np
import nrrd
from PIL.ImageOps import scale
from skimage.morphology import skeletonize_3d
from scipy.ndimage import label
import networkx as nx
import os
from scipy.ndimage import distance_transform_edt
from math import pi
import statistics
from scipy import stats
import pandas as pd


def avgBranchRadius(combinedGraph, radii):  # new_img_dst
    """Input: Built graph, distance transform of vesselmask
    First calculates additional distances of ending nodes from its background and adds it to the
    branch length. Next, for longer branches with distance > 5, pick intermediate pixels and calculate
    radial distance. Add length, radius, cross-sectional area to graph dataset. Repeat for smaller edges <= 5,
    only difference is that  intermediate pixels are only 3 pixels of the middle ones."""

    # get all ending nodes
    for u, v, o in combinedGraph.edges(data="object"):
        print(combinedGraph.edges())
        print(o)
        print(type(o))
        combinedGraph[u][v]["weight"] = len(o.pixels)
        print("#######branch radius estimation#########")

        combinedGraph[u][v]["type"] = "branchType"

        allpoints = []
        weight = combinedGraph[u][v]["weight"]
        allpoints.append(u)
        allpoints.append(o.pixels)
        allpoints.append(v)
        pixelDist = pixelDistance(allpoints)

        if weight >= 5:
            betweenPxlList = pointPicking(o.pixels, weight)
            csaAvg, radiusAvg, diameter = radialDistance(betweenPxlList, radii)  # new_img_dst
            combinedGraph[u][v][
                "length"] = pixelDist  # + endNodeAddition In my small vascular tree, I do not need endnodeadddition since the center line in my example start from center of first circle
            combinedGraph[u][v]["radius"] = radiusAvg
            combinedGraph[u][v]["diameter"] = diameter
            combinedGraph[u][v]["RatioLendiameter"] = pixelDist / diameter
            combinedGraph[u][v]["csa"] = csaAvg

        elif weight < 5 and weight >= 1:
            mid = weight // 2
            csaSmall, radiusSmall, diameter = radialDistance(o.pixels[mid - 1:mid + 1], radii)  # new_img_dst
            combinedGraph[u][v]["length"] = pixelDist  # + endNodeAddition
            combinedGraph[u][v]["radius"] = radiusSmall
            combinedGraph[u][v]["diameter"] = diameter
            combinedGraph[u][v]["RatioLendiameter"] = pixelDist / diameter
            combinedGraph[u][v]["csa"] = csaSmall

        else:

            combinedGraph[u][v]["length"] = 0
            combinedGraph[u][v]["radius"] = 0
            combinedGraph[u][v]["diameter"] = 0
            combinedGraph[u][v]["RatioLendiameter"] = 0
            combinedGraph[u][v]["csa"] = 0
            print("0 branch weight", weight)

            continue

    return combinedGraph


def radialDistance(middlePixels, radii):  # new_img_dst
    """Uses middle pixels and distance transform of vessel mask to return distance of centerline
    pixel from the nearest background, thereby calculating the radial distance of each vessel segment.
    Returns cross-sectional area and average radius from top 20% of vessel radius"""
    radiii = []
    for i in range(0, len(middlePixels)):
        R = middlePixels[i]
        radiii.append(radii[R])

    if len(radiii) != 0:
        radiii.sort(reverse=True)
        # toptrim = stats.trim1(radii,proportiontocut=0.1,tail='right')
        trimmedradius = stats.trim1(radiii, proportiontocut=0.9, tail='left')
        avg_radius = statistics.mean(trimmedradius)
        diameter = avg_radius * 2
    else:
        avg_radius = 0
        diameter = 0
    crossSecArea = pi * avg_radius ** 2
    return crossSecArea, avg_radius, diameter


def pixelDistance(pixelList):
    """From list of pixels, calculates euclidean distances between each pixels and returns total distance"""
    distance = 0
    a = pixelList
    allpairs = []
    allpairs.append(a[0])
    for i in range(len(a[1])):
        allpairs.append(a[1][i])
    allpairs.append(a[-1])
    pairs = [(allpairs[i], allpairs[i + 1]) for i in range(len(allpairs) - 1)]
    for p, q in pairs:
        p = pd.eval(p)
        q = pd.eval(q)
        distance = distance + np.linalg.norm(np.array(p) - np.array(q))
    return distance


def pointPicking(pixelList, weight):
    """returns list of pixels between 10-90% of whole list"""
    newList = pixelList[int(len(pixelList) * .0): int(len(pixelList) * .99)]
    return newList


def writeGraph(G, output_filepath):
    """Write edges of graph in xml format"""

    nx.write_edgelist(G, output_filepath, delimiter=',',
                      data=['radius', 'diameter', 'length', 'RatioLendiameter', 'ratioDiam', 'angle'])


def extract_centerline_points_with_connectivity(input_nrrd_file):
    # Read the NRRD file
    data, header = nrrd.read(input_nrrd_file)
    origin = np.array(header.get('space origin', [0.0, 0.0, 0.0]))
    spacing = np.array([header.get(f'space directions {i}', 1.0) for i in range(3)])

    # Skeletonize the vessel mask
    skeleton = skeletonize_3d(data > 0)

    # Label connected components of the skeleton
    labeled_skeleton, num_features = label(skeleton)

    # Create a graph from the skeleton
    G = nx.Graph()
    shape = skeleton.shape
    for x in range(shape[0]):
        for y in range(shape[1]):
            for z in range(shape[2]):
                if labeled_skeleton[x, y, z] > 0:
                    G.add_node((x, y, z))
                    for dx, dy, dz in [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]:
                        neighbor_x, neighbor_y, neighbor_z = x + dx, y + dy, z + dz
                        if 0 <= neighbor_x < shape[0] and 0 <= neighbor_y < shape[1] and 0 <= neighbor_z < shape[2]:
                            if labeled_skeleton[neighbor_x, neighbor_y, neighbor_z] > 0:
                                G.add_edge((x, y, z), (neighbor_x, neighbor_y, neighbor_z))

    # Compute the distance transform of the vessel mask
    distance_transform = distance_transform_edt(data > 0) * spacing.min()

    # Convert pixel coordinates to physical coordinates
    centerline_points = [
        tuple(float(coord) for coord in origin + np.multiply(node, spacing)) for node in G.nodes
    ]
    connectivity = [
        (
            tuple(float(coord) for coord in origin + np.multiply(edge[0], spacing)),
            tuple(float(coord) for coord in origin + np.multiply(edge[1], spacing))
        )
        for edge in G.edges
    ]

    # Calculate radius for each edge
    radii = []
    for edge in G.edges:
        node1, node2 = edge
        midpoint = [(n1 + n2) / 2 for n1, n2 in zip(node1, node2)]  # Midpoint in pixel space
        radius = distance_transform[int(midpoint[0]), int(midpoint[1]), int(midpoint[2])]
        radii.append(radius)

    scale = np.array(radii)
    sigma0 = 1 / np.sqrt(2.) / 2.
    selfspacing = np.prod(spacing) ** (1 / 3.0)
    sigmap = 1 / np.sqrt(2.) / 2.
    mask = scale < (2. / np.sqrt(2) * sigma0)
    rad = np.zeros(mask.shape)
    rad[mask] = np.sqrt(2.) * (np.sqrt((scale[mask] * selfspacing) ** 2.0 + (
            sigma0 * selfspacing) ** 2.0) - 0.5 * sigma0 * selfspacing)
    rad[~mask] = np.sqrt(2.) * (np.sqrt((scale[~mask] * selfspacing) ** 2.0 + (
            sigmap * selfspacing) ** 2.0) - 0.5 * sigmap * selfspacing)

    new_dic = {key: rad[i] for i, key in enumerate(radii)}
    full_graph = avgBranchRadius(G, new_dic)

    return centerline_points, connectivity, radii, rad, full_graph


# Example usage
input_nrrd = os.path.join("C:\\Users\GGPC\\Downloads\\P-7B5VH\\P-7B5VH\\Pre\\Vessel\\MPA_mask\\Segmentation.nrrd")
centerline_points, connectivity, radii, rad, full_graph = extract_centerline_points_with_connectivity(input_nrrd)
writeGraph(full_graph, os.path.join("C:\\Users\GGPC\\Downloads\\P-7B5VH\\P-7B5VH\\Pre\\Vessel\\MPA_mask\\TEST.csv"))

# Save or process the results
print("Centerline Points:", centerline_points)
print("Connectivity:", connectivity)
print("Radii:", radii)
