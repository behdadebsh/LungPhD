import numpy as np

map_file_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/RML_mapping.txt'
mesh_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Vessel/CTEPH10_Artery_Full'
# intensity_path = '/hpc/bsha219/lung/Data/CTEPH/Arthur/INTENSITY-BASED_FLOW_MAPS/proximal_flow_mapping/CTEPH3/CTEPH3_flow_diff_fractions_avg_0g.exelem'
cluster_file = '/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/RML_cluster_map_terminals.exnode'
grid_file = '/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/RML_DataGrid.ipdata'
# num_terminal = get_num_terminals(map_file_path)
# terminals = get_terminal_array(map_file_path, num_terminal)


def search_string_in_file(file_name, string_to_search):
    """Search for the given string in file and return lines containing that string,
    along with line numbers"""
    line_number = 0
    list_of_results = []
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            line_number += 1
            if string_to_search in line and 'Element' not in line:
                # If yes, then add the line number & line as a tuple in the list
                list_of_results.append((line_number, line.rstrip()))
    # Return list of tuples containing line numbers and lines where string is found
    return list_of_results

with open(map_file_path) as f_in:
    lines = (line.rstrip() for line in f_in)
    map_contents = list(line for line in lines if line)
num_terminals = len(map_contents)

with open(mesh_path + '.ipelem') as f_in:
    lines = (line.rstrip() for line in f_in)
    elem_contents = list(line for line in lines if line)

with open(cluster_file) as f_in:
    lines = (line.rstrip() for line in f_in)
    cluster_contents = list(line for line in lines if line)

with open(grid_file) as f_in:
    lines = (line.rstrip() for line in f_in)
    grid_contents = list(line for line in lines if line)

term_field = np.zeros((num_terminals, 6), dtype=float) # includes data point, terminal number, x,y,z and cluster number

for line in range(num_terminals):
    term_field[line][0] = map_contents[line].split()[0] # data point number extracted
    elem = int(map_contents[line].split()[1]) # find the corresponding elem number to the data point
    counter = 0
    for elem_line in range(len(elem_contents)):
        if elem_contents[elem_line].split()[0] == 'Element' and int(elem_contents[elem_line].split()[-1]) == elem:
            # print("here!")
            elem_term_node = elem_contents[elem_line + 5].split()[-2]
            matched = search_string_in_file(mesh_path + '.ipelem', ' ' + str(elem_term_node))
            for rep in matched:
                if str(elem_term_node) in rep[1].split():
                    counter = counter + 1
            if counter > 1:
                elem_term_node = elem_contents[elem_line + 5].split()[-1]
            term_field[line][1] = elem_term_node
    for cluster_line in range(len(cluster_contents)):
        if(cluster_contents[cluster_line].split()[0]=='Node:' and cluster_contents[cluster_line].split()[-1] == elem_term_node):
            # term_field[line][2] = cluster_contents[cluster_line + 1].split()[-1]
            # term_field[line][3] = cluster_contents[cluster_line + 2].split()[-1]
            # term_field[line][4] = cluster_contents[cluster_line + 3].split()[-1]
            term_field[line][5] = cluster_contents[cluster_line + 5].split()[-1]
            # print(term_field[line][5])
            # print(type(term_field[line][5]))
    for grid_line in range(1, len(grid_contents)):
        # print(float(grid_contents[grid_line].split()[0]))
        # print(term_field[line][0])
        if (float(grid_contents[grid_line].split()[0]) - 10000 == term_field[line][0]):  ############################# ATTENTION!!!!!! #############################################
            # print('here')
            term_field[line][2] = grid_contents[grid_line].split()[1]
            term_field[line][3] = grid_contents[grid_line].split()[2]
            term_field[line][4] = grid_contents[grid_line].split()[3]


with open("/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/RML_CL1.ipdata", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
    for i in range(len(term_field)):
        if (term_field[i][5] == 7):
            f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
                                                  float(term_field[i][3]), float(term_field[i][4])))

with open("/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/RML_CL2.ipdata", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
    for i in range(len(term_field)):
        if (term_field[i][5] == 8):
            f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
                                                  float(term_field[i][3]), float(term_field[i][4])))

with open("/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/RML_CL3.ipdata", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
    for i in range(len(term_field)):
        if (term_field[i][5] == 9):
            f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
                                                  float(term_field[i][3]), float(term_field[i][4])))

with open("/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/RML_CL4.ipdata", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
    for i in range(len(term_field)):
        if (term_field[i][5] == 10):
            f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
                                                  float(term_field[i][3]), float(term_field[i][4])))

# with open("/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/RLL_CL5.ipdata", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
#     for i in range(len(term_field)):
#         if (term_field[i][5] == 15):
#             f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
#                                                   float(term_field[i][3]), float(term_field[i][4])))
#
# with open("/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/RLL_CL6.ipdata", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
#     for i in range(len(term_field)):
#         if (term_field[i][5] == 16):
#             f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
#                                                   float(term_field[i][3]), float(term_field[i][4])))

# with open("/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/LUL_CL7.ipdata", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
#     for i in range(len(term_field)):
#         if (term_field[i][5] == 29):
#             f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
#                                                   float(term_field[i][3]), float(term_field[i][4])))
#
# with open("/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/LUL_CL8.ipdata", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
#     for i in range(len(term_field)):
#         if (term_field[i][5] == 30):
#             f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
#                                                   float(term_field[i][3]), float(term_field[i][4])))




