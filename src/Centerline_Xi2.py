#!/usr/bin/env python

import numpy as np
from opencmiss.zinc.context import Context
from opencmiss.zinc.status import OK as ZINC_OK
from utils import writeExDataFile

# This script is written to calculate a centerline points (including the junction point) and radii respective to those
# points at each point exploiting the *.exnode and .*exelem files reading them using OpenCMISS.Zinc and respective to
# Xi_2 direction which in this case is towards the length of the tubes for a single bifurcation template. The centerline
#  points are written as *.exdata file format which is compatible with CMGUI.

def getCoordinatesForXi2(el_info, el_ids, xi_2_index):
    """
    Get the coordinates for the ring specified by the el_ids at the given xi_2 location using the
    supplied el_info.
    :param el_info: Dictionary of el_ids and a list of coordinates increasing in xi_1 fastest.
    :param el_ids: A list of element ids that describe the ring.
    :param xi_2_index: The index of the xi_2 location to extract the coordinates for.
    :return: A list of coordinates that form the ring.
    """
    coordinate_list = []
    for el_id in el_ids:
        all_coordinates_for_element = el_info[el_id]
        coordinate_list.extend(all_coordinates_for_element[5*xi_2_index:5*(xi_2_index + 1)])
    return coordinate_list


# def junctionCoordinate(element_locations, ring_numbers):
#     """
#     Get the coordinates for the junction of the bifurcation using
#     supplied el_info.
#     :param element_locations: Dictionary of el_ids and a list of coordinates increasing in xi_1 fastest.
#     :param rings: A list of element ids that describe the ring.
#     :return: A list of coordinates that form the ring.
#     """
#     summation = getCoordinatesForXi2(element_locations, ring_numbers[0], 4) + getCoordinatesForXi2(element_locations, ring_numbers[1], 0) + getCoordinatesForXi2(element_locations, ring_numbers[2], 0)
#     junction_point = np.mean(summation)
#     return junction_point


context = Context("Centerline")
region = context.getDefaultRegion()
directory = '/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/SurfaceFEMesh/'
status = region.readFile(directory + "MPA_fitted.exnode")
status = region.readFile(directory + "MPA_fitted.exelem")
fieldmodule = region.getFieldmodule()
field = fieldmodule.findFieldByName("coordinates")
cache = fieldmodule.createFieldcache()

center_points = []
radii = []
ring_1 = [1, 2, 3, 4]
ring_2 = [25, 26, 27, 28]
ring_3 = [5, 6, 7, 8]
ring_4 = [9, 10, 11, 12]
ring_5 = [13, 14, 15, 16]
ring_6 = [17, 18, 19, 20]
ring = [ring_1, ring_2, ring_3, ring_4, ring_5, ring_6]
for i in range(6):
    if i == 4 or i == 5:
        xi1 = [0, 0.2, 0.4, 0.6, 0.8]
        xi2 = [0, 0.25, 0.5, 0.75, 1]
    else:
        xi1 = [0, 0.2, 0.4, 0.6, 0.8]
        xi2 = [0, 0.25, 0.5, 0.75]
    mesh = fieldmodule.findMeshByDimension(2)
    el_iter = mesh.createElementiterator()
    element = el_iter.next()

    el_locations = {}
    while element.isValid():
        el_id = element.getIdentifier()
        el_locations[el_id] = []
        for xi_2 in xi2:
            for xi_1 in xi1:
                xi = [xi_1, xi_2]
                cache.setMeshLocation(element, xi)
                result, outValues = field.evaluateReal(cache, 3)
                # Check result for errors, Use outValues
                if result == ZINC_OK:
                    # print(element.getIdentifier(), outValues)
                    el_locations[el_id].append(outValues)
                else:
                    break
        element = el_iter.next()
    for xi_2 in xi2:
        ring_coordinates = getCoordinatesForXi2(el_locations, ring[i], xi2.index(xi_2))
        # print('ring is:', ring_coordinates)
        center_point = np.mean(ring_coordinates, axis=0)
        radius = []
        for nodes_on_ring in range(len(ring_coordinates)):
            r = np.sqrt((ring_coordinates[nodes_on_ring][0] - center_point[0])**2 + (ring_coordinates[nodes_on_ring][1] - center_point[1])**2 + (ring_coordinates[nodes_on_ring][2] - center_point[2])**2)
            radius.append(r)
        radii.append(np.mean(radius))
            # radii.append(radii)
        # print('centerx is:', center_point[0])
        # print(center_point)
        center_points.append(center_point)
c = getCoordinatesForXi2(el_locations, ring_2, 4)
d = getCoordinatesForXi2(el_locations, ring_3, 0)
e = getCoordinatesForXi2(el_locations, ring_4, 0)
junction_rings = c + d + e
junction_point = np.mean(junction_rings, axis=0)
center_points.append(junction_point)
junction_radii = []
for j in range(len(junction_rings)):
    junction_radius = np.sqrt((junction_rings[j][0] - junction_point[0]) ** 2 + (junction_rings[j][1] - junction_point[1]) ** 2 + (junction_rings[j][2] - junction_point[2]) ** 2)
    junction_radii.append(junction_radius)
junction_radii = np.mean(junction_radii)
radii.append(junction_radii)
writeExDataFile('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/1DMesh/centre_points.exdata', center_points, radii)
