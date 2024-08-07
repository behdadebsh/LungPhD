import math
import numpy as np
import os, sys


def read_nodes(name):

    # read a basic ipnode file (no derivatives)

    nodefile = open(name + '.ipnode', 'r')
    node_coordinates = {}
    line = nodefile.readline()
    while line:
        if 'node number' in line.lower():
            node = int(line.split(']:')[1])
            line = nodefile.readline()
            x_coord = float(line.split(']:')[1])
            line = nodefile.readline()
            y_coord = float(line.split(']:')[1])
            line = nodefile.readline()
            z_coord = float(line.split(']:')[1])
            node_coordinates[node] = np.array([x_coord, y_coord, z_coord])
        line = nodefile.readline()

    nodefile.close()

    return node_coordinates


def read_elements(name):
    # read a basic 1d element file (no derivatives)

    elemfile = open(name + '.ipelem', 'r')
    element_nodes = {}
    line = elemfile.readline()
    while line:
        if 'Element number' in line:
            elem = int(line.split()[4])
            for i in range(5): line = elemfile.readline()
            node1 = int(line.split()[8])
            node2 = int(line.split()[9])
            element_nodes[elem] = [node1, node2]
        line = elemfile.readline()

    elemfile.close()

    return element_nodes


def read_fields(name):

    # read basic field file, radii
    fieldfile = open(name + '.ipfiel', 'r')
    element_field = {}
    line = fieldfile.readline()
    while line:
        if 'Element number' in line:
            elem = int(line.split()[2])
            line = fieldfile.readline()
            radius = float(line.split(']:')[1])
            element_field[elem] = radius
        line = fieldfile.readline()

    fieldfile.close()

    return element_field


def connect_tree(element_nodes):
    # record which elements are connected to which
    connect_down = {}
    connect_up = {}
    for elem1 in element_nodes:
        inlet_node1 = element_nodes[elem1][0]
        outlet_node1 = element_nodes[elem1][1]
        downstream = []
        for elem2 in element_nodes:
            inlet_node2 = element_nodes[elem2][0]
            outlet_node2 = element_nodes[elem2][1]
            if outlet_node1 == inlet_node2:
                downstream.append(elem2)
            if outlet_node2 == inlet_node1:
                connect_up[elem1] = elem2
        if downstream:  # not empty list
            connect_down[elem1] = downstream

    return connect_up, connect_down


def largest_inlet(connect_up, connect_down, radii):
    # get the inlet or outlet vessel with largest radius value
    radius_largest = -1.0e+6
    for elem in radii:
        if elem not in connect_up:
            # this is an inlet
            if radii[elem] > radius_largest:
                elem_largest = elem
                radius_largest = radii[elem]
        if elem not in connect_down:
            # this is an outlet; also check just in case
            if radii[elem] > radius_largest:
                elem_largest = elem
                radius_largest = radii[elem]

    return elem_largest


def airway_list(target_key, connect_down):
    elem_list = []
    elem_list.append(target_key)
    elem_count = len(connect_down.get(target_key, []))  # number of elements in the list
    while elem_count == 1:  # only attached to one element. refined branch
        target_key = connect_down[target_key][0]  # get the next element number
        elem_list.append(target_key)
        elem_count = len(connect_down.get(target_key, []))  # number of elements in the list

    return elem_list


def calculate_angle(vect1, vect2):

    # angle (in radians) between two given (non-unit) vectors

    angle_rads = np.arccos(np.dot(vect1, vect2) / (np.linalg.norm(vect1)*np.linalg.norm(vect2)))
    return angle_rads


def child_labels(parent, child1, child2, airways, connect_down, element_nodes, coords, cf_vect):

    elem0 = airways[parent][-1] # last element in list
    elem1 = connect_down[elem0][0]
    elem2 = connect_down[elem0][1]
    node0 = element_nodes[elem1][0]
    node1 = element_nodes[elem1][1]
    node2 = element_nodes[elem2][1]
    vect1 = coords[node1] - coords[node0]
    vect2 = coords[node2] - coords[node0]
    angle1 = calculate_angle(vect1, cf_vect)
    angle2 = calculate_angle(vect2, cf_vect)
    if angle1 < angle2:
        airways[child1] = airway_list( elem1, connect_down )
        airways[child2] = airway_list( elem2, connect_down )
    else:
        airways[child2] = airway_list( elem1, connect_down )
        airways[child1] = airway_list( elem2, connect_down )

    return airways[child1], airways[child2]


def label_arteries(MPA, element_nodes, connect_down, coords):

    x_vect = np.array([1,0,0])
    y_vect = np.array([0,1,0])
    z_vect = np.array([0,0,1])
    arteries = {}
    arteries['MPA'] = airway_list( MPA, connect_down)

    arteries['RMA'], arteries['LMA'] = child_labels('MPA', 'RMA', 'LMA', arteries, connect_down, element_nodes, coords, x_vect )
    arteries['RLMA'], arteries['RULA'] = child_labels('RMA', 'RLMA', 'RULA', arteries, connect_down, element_nodes, coords, -z_vect )
    arteries['RMLA'], arteries['RLLA'] = child_labels('RLMA', 'RMLA', 'RLLA', arteries, connect_down, element_nodes, coords, y_vect )
    arteries['LLMA'], arteries['LULA'] = child_labels('LMA', 'LLMA', 'LULA', arteries, connect_down, element_nodes, coords, -z_vect )
    arteries['ling'], arteries['LLLA'] = child_labels('LLMA', 'ling', 'LLLA', arteries, connect_down, element_nodes, coords, y_vect )

    return arteries


def annotate_arteries(upper_nodes, upper_elems, upper_radii):
    node_coordinates = read_nodes(upper_nodes)
    element_nodes = read_elements(upper_elems)
    radii = read_fields(upper_radii)
    connect_up, connect_down = connect_tree(element_nodes)
    MPA = largest_inlet(connect_up, connect_down, radii)
    arteries = label_arteries(MPA, element_nodes, connect_down, node_coordinates)

    return arteries