import numpy as np
import matplotlib.pyplot as plt
# Though the following import is not directly being used, it is required
# for 3D projection to work
from mpl_toolkits.mplot3d import Axes3D
import time
from sklearn.cluster import KMeans


def match_terminal(map_file_path, mesh_path, cluster_file_path):
    """
    This function uses the ipnode and ipelem files to match cluster labels and datagrid points within a lobe
    :param map_file_path: data to terminal element map file - output from growing
    :param mesh_path: path to ipnode and ipelem file
    :param cluster_file: path to terminal exnode file with cluster and intensity fields (RLL_cluster_map_terminals.exnode)
    :return: write_to_file data point and its corresponding cluster
    """
    with open(map_file_path) as f_in:
        lines = (line.rstrip() for line in f_in)
        map_contents = list(line for line in lines if line)
    num_terminals = len(map_contents) - 1
    print(num_terminals)
    # with open(mesh_path + '.ipelem') as f_in:
    #     lines = (line.rstrip() for line in f_in)
    #     elem_contents = list(line for line in lines if line)
    # with open(cluster_file_path + '.exnode') as f_in:
    #     lines = (line.rstrip() for line in f_in)
    #     node_contents = list(line for line in lines if line)
    # num_terminals = get_num_terminals(map_file_path)
    # term_field = np.zeros((num_terminals, 4), dtype=float) # includes x,y,z and cluster number
    # term_nodes = np.zeros(num_terminals)
    # term_elems = get_terminal_array(map_file_path,num_terminals)
    #
    # i = 0  # num_terminal counter
    # for elem in term_elems:
    #     counter = 0
    #     for elem_line in range(len(elem_contents)):
    #         if elem_contents[elem_line].split()[0] == 'Element' and float(elem_contents[elem_line].split()[-1]) == elem:
    #             elem_term_node = elem_contents[elem_line+5].split()[-2]
    #             matched = search_string_in_file(mesh_path + '.ipelem', ' ' + str(elem_term_node))
    #             for rep in matched:
    #                 if str(elem_term_node) in rep[1].split():
    #                     counter = counter + 1
    #             if counter > 1:
    #                 elem_term_node = elem_contents[elem_line + 5].split()[-1]
    #             term_nodes[i] = elem_term_node
    #     for node_line in range(len(node_contents)):
    #         if node_contents[node_line].split()[0] == 'Node' and node_contents[node_line].split()[-1] == elem_term_node:
    #             if node_contents[node_line+2].split()[0] == 'For':
    #                 x = node_contents[node_line+3].split()[-1]
    #                 y = node_contents[node_line+8].split()[-1]
    #                 z = node_contents[node_line+13].split()[-1]
    #                 term_field[i][0] = float(x)
    #                 term_field[i][2] = float(y)
    #                 term_field[i][3] = float(z)
    #             else:
    #                 x = node_contents[node_line+2].split()[-1]
    #                 y = node_contents[node_line+4].split()[-1]
    #                 z = node_contents[node_line+6].split()[-1]
    #                 term_field[i][0] = float(x)
    #                 term_field[i][2] = float(y)
    #                 term_field[i][3] = float(z)



def read_flow_intensity_map(map_file_path, mesh_path, num_terminal, terminals):
    """
    read flow intensity mappings from file on terminal elements and x,y,z coordinates of terminal node of respective
    element stored a numpy array
    :param map_file_path: path to flow diff fraction file
    :param mesh_path: path to node and elem file
    :param num_terminal: number of terminals
    :param terminals: array of terminal elements
    :return: term_field: an array consisting of x,y,z and flow_diff_mapping in respective order
    """
    term_field = np.zeros((num_terminal, 4), dtype=float)
    term_nodes = np.zeros(num_terminal)

    with open(map_file_path) as f_in:
        lines = (line.rstrip() for line in f_in)
        map_contents = list(line for line in lines if line)
    with open(mesh_path + '.ipelem') as f_in:
        lines = (line.rstrip() for line in f_in)
        elem_contents = list(line for line in lines if line)
    with open(mesh_path + '.ipnode') as f_in:
        lines = (line.rstrip() for line in f_in)
        node_contents = list(line for line in lines if line)
        # intensity_contents = list(line for line in lines if line)
    i = 0  # num_terminal counter
    for elem in terminals:
        counter = 0
        for elem_line in range(len(elem_contents)):
            if elem_contents[elem_line].split()[0] == 'Element' and float(elem_contents[elem_line].split()[-1]) == elem:
                elem_term_node = elem_contents[elem_line+5].split()[-2]
                matched = search_string_in_file(mesh_path + '.ipelem', ' ' + str(elem_term_node))
                for rep in matched:
                    if str(elem_term_node) in rep[1].split():
                        counter = counter + 1
                if counter > 1:
                    elem_term_node = elem_contents[elem_line + 5].split()[-1]
                term_nodes[i] = elem_term_node
        for node_line in range(len(node_contents)):
            if node_contents[node_line].split()[0] == 'Node' and node_contents[node_line].split()[-1] == elem_term_node:
                if node_contents[node_line+2].split()[0] == 'For':
                    x = node_contents[node_line+3].split()[-1]
                    y = node_contents[node_line+8].split()[-1]
                    z = node_contents[node_line+13].split()[-1]
                    term_field[i][0] = float(x)
                    term_field[i][2] = float(y)
                    term_field[i][3] = float(z)
                else:
                    x = node_contents[node_line+2].split()[-1]
                    y = node_contents[node_line+4].split()[-1]
                    z = node_contents[node_line+6].split()[-1]
                    term_field[i][0] = float(x)
                    term_field[i][2] = float(y)
                    term_field[i][3] = float(z)
        for map_line in range(len(map_contents)):
            if map_contents[map_line].split()[0] == 'Element:':
                if float(map_contents[map_line].split()[1]) == elem:
                    term_field[i][1] = float(map_contents[map_line+2].split()[0])
        i = i + 1
    return term_field, term_nodes


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


def get_num_terminals(file_path):
    """
    Get the number of terminals for the 1D mesh of the subject based on map.
    This function will only get the number of terminals after growing without the need of running simulations.
    :param file_path: path to data to terminal map
    :return: Number of terminals
    """
    map_file = open(file_path)
    file_content = map_file.readlines()
    num_terminals = len(file_content)
    return num_terminals

def get_terminal_array(file_path, num_terminal):
    """
    Get the array of terminal element numbers
    :param file_path: Path to data to terminal map
    :param num_terminal: number of terminals for array size allocation
    :return: terminals: an array of terminal elems
    """
    terminals = np.zeros(num_terminal)
    map_file = open(file_path)
    file_content = map_file.readlines()
    for line in range(len(file_content)):
        terminals[line] = file_content[line].split()[1]
    return terminals


def main():
    # start = time.clock()
    # map_file_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH1/FRC/Lobe/RML_data_node_map.out'
    # map_file_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH4/FRC/Lobe/CTEPH4_whole_map.txt'
    # map_file_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/LLL_mapping.txt'
    map_file_path = '/hpc/bsha219/Research related/EMBC2023/Intensity_mapping/mapping.txt'
    # mesh_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Vessel/CTEPH10_Artery_Full'
    mesh_path = '/hpc/bsha219/Research related/EMBC2023/Geometry/block_grown'
    # intensity_path = '/hpc/bsha219/lung/Data/CTEPH/Arthur/INTENSITY-BASED_FLOW_MAPS/proximal_flow_mapping/CTEPH3/CTEPH3_flow_diff_fractions_avg_0g.exelem'
    # intensity_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Intensity_mapping/CTEPH10_flow_diff_fractions_avg_RM1.exelem'
    intensity_path = '/hpc/bsha219/Research related/EMBC2023/Intensity_mapping/flow_diff_fraction.exelem'
    num_terminal = get_num_terminals(map_file_path)
    terminals = get_terminal_array(map_file_path, num_terminal)
    X, terminal_nodes = read_flow_intensity_map(intensity_path, mesh_path, num_terminal, terminals)

    sacre_bleu = KMeans(n_clusters=16).fit(X)
    labels = sacre_bleu.labels_
    print(labels)
    clusters = sacre_bleu.cluster_centers_
    print(clusters)
    fig = plt.figure(1, figsize=(4, 3))
    ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
    ax.scatter(X[:, 3], X[:, 0], X[:, 2], c=labels.astype(float), edgecolor='k')
    # end = time.clock()
    # print("elapsed time:", end - start)

    # find the underperfused clusters
    for cluster in range(len(clusters)):
        if clusters[cluster][1] < 0:  # underperfused
            print("Underperfused:", cluster+1, "Avg. Intensity:", clusters[cluster][1])

    termdict = []

    for i in range(len(terminal_nodes)):
        termdict.append(dict(node_num=terminal_nodes[i], x=X[i][0], y=X[i][2], z=X[i][3], intensity=X[i][1],
                             cluster=labels[i]+1))

    termdict = sorted(termdict, key=lambda k: k['node_num'])

    # Lobe = 'whole20'
    Lobe = 'block'
    print(Lobe)
    filename = Lobe + '_cluster_map_terminals.exnode'
    export_path = '/hpc/bsha219/Research related/EMBC2023/Intensity_mapping/'
    with open(export_path + filename, 'w') as f:
        f.write(' Group name : terminal clusters\n')
        f.write(' #Fields=3\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('  x.  Value index=1, #Derivatives=0\n')
        f.write('  y.  Value index=2, #Derivatives=0\n')
        f.write('  z.  Value index=3, #Derivatives=0\n')
        f.write(' 2) intensity, field, rectangular cartesian, #Components=1\n')
        f.write('  1.  Value index=4, #Derivatives=0\n')
        f.write(' 3) cluster, field, rectangular cartesian, #Components=1\n')
        f.write('  1.  Value index=5, #Derivatives=0\n')

        for i in range(len(X)):
            f.write(' Node:           %d\n' % (termdict[i]['node_num']))
            f.write('     %s\n' % float(termdict[i]['x']))
            f.write('     %s\n' % float(termdict[i]['y']))
            f.write('     %s\n' % float(termdict[i]['z']))
            f.write('     %s\n' % (termdict[i]['intensity']))
            f.write('     %s\n' % (termdict[i]['cluster']))

    plt.show()


if __name__ == '__main__':
    main()
