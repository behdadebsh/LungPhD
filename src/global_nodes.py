import os
import numpy as np
import re

def global_nodes(Vessel_Path, Artery_File, Number_Elements):
    # ELEMENT GLOBAL NODES
    # the file CTEPH4_Artery_Full.exelem contains the initial arterial mesh
    # with just element to node mapping data

    os.chdir(Vessel_Path)
    nodes = np.zeros((Number_Elements, 2))
    with open(Artery_File, 'r') as fid:
        content = fid.read()

    # Extract all numbers after 'Nodes:'
    matches = re.findall(r'Nodes:\s*(\d+)\s+(\d+)', content)

    # Convert the matches to a 2D NumPy array
    nodes = np.array([list(map(int, match)) for match in matches])
    elem_global_nodes = nodes

    # PARENT CHILD ARRAY
    l = len(elem_global_nodes)
    max_index = np.max(elem_global_nodes)
    node_par = np.zeros((max_index+1, 1), dtype=int)    # global nodes as parent/child markers
    elem_children = np.zeros((l, 2), dtype=int)
    node_par[0, 0] = 0    # global node 1 is inlet so N/A
    node_par[1, 0] = 1    # global node 2 is local node 2 of element 1 -> element 1 is stored as parent of this node

    for k in range(1, l):
        if elem_children[node_par[elem_global_nodes[k, 0], 0], 0] == 0:
            elem_children[node_par[elem_global_nodes[k, 0], 0], 0] = k+1
        else:
            elem_children[node_par[elem_global_nodes[k, 0], 0], 1] = k+1
        node_par[elem_global_nodes[k, 1], 0] = k

    # ELEMENT ORDER FROM GLOBAL NODES
    # This is just to specify an order for the calculation of impedance and
    # pressure so not a big deal that each element in the MPA has its own order

    max_index = np.max(elem_global_nodes)
    node_gen = np.zeros((max_index+1, 1), dtype=int)    # array for generation signaled by local node 1
    elem_gen = np.zeros((l, 1), dtype=int)
    node_gen[0, 0] = 1    # not used but for bookkeeping
    node_gen[1, 0] = 2    # if element local node 1 = node 2 => gen = 2
    elem_gen[0, 0] = 1    # set inlet element as gen 1

    for k in range(1, l):
        elem_gen[k, 0] = node_gen[elem_global_nodes[k, 0], 0]+2
        node_gen[elem_global_nodes[k, 1], 0] = node_gen[elem_global_nodes[k, 0], 0]+1

    return elem_global_nodes, elem_gen, elem_children