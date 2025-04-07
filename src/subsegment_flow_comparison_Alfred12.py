import os
import numpy as np


def get_num_elements(artery_file):
    with open(artery_file) as f_in:
        lines = (line.rstrip() for line in f_in)
        file_content = list(line for line in lines if line)
    for file_line in (file_content):
        if (file_line.split()[0]) == 'Element:':
            num_elems = int(file_line.split()[1])
    return num_elems


def adjust_hu_values(node_list):
    # Extract all HU values from the dictionary
    min_hu = min(node['HU'] for node in node_list)
    print(min_hu)

    # Add the minimum HU value to all nodes
    for node in node_list:
        node['HU'] += min_hu

    return node_list


root_folder = '/hpc/bsha219/lung/Data'
study = 'CTEPH'
protocol = 'Pre'
subject_name = 'Alfred12'
subject_path = os.path.join(root_folder, study, subject_name, protocol)
flow_file = os.path.join(subject_path, 'Perfusion', 'Baseline0.0') + '/' + subject_name + '_flow_perf.exelem'
num_elems = get_num_elements(flow_file)
flows = np.zeros((num_elems, 1))  # preallocate array
with open(flow_file) as f_in:
    lines = (line.rstrip() for line in f_in)
    flow_file_content = list(line for line in lines if line)
i = 0
for line in range(len(flow_file_content)):
    if flow_file_content[line].split()[0] == 'Element:' and flow_file_content[line].split()[1] == '1':
        inlet_flow = float(flow_file_content[line+2].split()[0])
        flows[i] = float(flow_file_content[line+2].split()[0])
        i = i+1
    elif flow_file_content[line].split()[0] == 'Element:':
        flows[i] = float(flow_file_content[line+2].split()[0])
        i = i+1

normalised_flows = flows/inlet_flow

# cluster mapping dictionary based on subsegments
branch_mapping = {
    1: {'Lobe': 'RUL', 'branch_name': 'Anterior_lower', 'element': 35},
    2: {'Lobe': 'RUL', 'branch_name': 'Apical_upper', 'element': 49},
    3: {'Lobe': 'RUL', 'branch_name': 'Apical_lower', 'element': 50},
    4: {'Lobe': 'RUL', 'branch_name': 'Posterior_lower', 'element': 52},
    5: {'Lobe': 'RUL', 'branch_name': 'Posterior_upper', 'element': 51},
    6: {'Lobe': 'RUL', 'branch_name': 'Anterior_upper', 'element': 36},
    7: {'Lobe': 'RML', 'branch_name': 'Medial_anterior', 'element': 57},
    8: {'Lobe': 'RML', 'branch_name': 'Lateral_anterior', 'element': 59},
    9: {'Lobe': 'RML', 'branch_name': 'Medial_posterior', 'element': 58},
    10: {'Lobe': 'RML', 'branch_name': 'Lateral_posterior', 'element': 60},
    11: {'Lobe': 'RLL', 'branch_name': 'Superior_upper', 'element': 54},
    12: {'Lobe': 'RLL', 'branch_name': 'Lateral', 'element': 65},
    13: {'Lobe': 'RLL', 'branch_name': 'Superior_lower', 'element': 53},
    14: {'Lobe': 'RLL', 'branch_name': 'Medial', 'element': 68},
    15: {'Lobe': 'RLL', 'branch_name': 'Posterior', 'element': 67},
    16: {'Lobe': 'RLL', 'branch_name': 'Anterior', 'element': 66},
    17: {'Lobe': 'LLL', 'branch_name': 'Anterior', 'element': 63},
    18: {'Lobe': 'LLL', 'branch_name': 'Lateral', 'element': 64},
    19: {'Lobe': 'LLL', 'branch_name': 'Posterior', 'element': 62},
    20: {'Lobe': 'LLL', 'branch_name': 'Superior_upper', 'element': 45},
    21: {'Lobe': 'LLL', 'branch_name': 'Medial', 'element': 61},
    22: {'Lobe': 'LLL', 'branch_name': 'Superior_lower', 'element': 46},
    23: {'Lobe': 'LUL', 'branch_name': 'Apical_lower', 'element': 41},
    24: {'Lobe': 'LUL', 'branch_name': 'Apical_upper', 'element': 42},
    25: {'Lobe': 'LUL', 'branch_name': 'Anterior_lower', 'element': 27},
    26: {'Lobe': 'LUL', 'branch_name': 'Lingular_lower', 'element': 30},
    27: {'Lobe': 'LUL', 'branch_name': 'Anterior_upper', 'element': 28},
    28: {'Lobe': 'LUL', 'branch_name': 'Posterior_lower', 'element': 43},
    29: {'Lobe': 'LUL', 'branch_name': 'Lingular_upper', 'element': 29},
    30: {'Lobe': 'LUL', 'branch_name': 'Posterior_upper', 'element': 44},
}
cluster_directory = os.path.join(subject_path, 'Intensity_mapping')
lobes = ['RUL', 'RML', 'RLL', 'LUL', 'LLL']
nodes = []
for lobe in lobes:
    with open(os.path.join(cluster_directory, lobe + '_cluster_map_terminals' + '.exnode')) as lobe_file:
        lines = (line.rstrip() for line in lobe_file)
        file_content = list(line for line in lines if line)
    for line in range(len(file_content)):
        if file_content[line].split()[0] == 'Node:':
            nodes.append({'x': float(file_content[line+1]), 'y': float(file_content[line+2]), 'z': float(file_content[line+3]),
                          'HU': float(file_content[line+4]), 'cluster': int(file_content[line+5])})

adjusted_node_dict = adjust_hu_values(nodes)
# Calculate the total HU value across all nodes
total_HU = sum(float(node['HU']) for node in adjusted_node_dict)
cluster_HU_sums = {}

for node in adjusted_node_dict:
    cluster = node['cluster']
    if cluster not in cluster_HU_sums:
        cluster_HU_sums[cluster] = 0
    cluster_HU_sums[cluster] += float(node['HU'])

normalised_cluster_HU = {cluster: hu_sum / total_HU for cluster, hu_sum in cluster_HU_sums.items()}
# compare flow
# sum = 0
output_file_path = os.path.join(subject_path, 'Intensity_mapping', 'flow_comparison_NEW.log')
with open(output_file_path, 'w') as fid:
    for cluster, normalised_HU in normalised_cluster_HU.items():
        if cluster in branch_mapping:
            element = branch_mapping[cluster]['element']
            # sum = sum + normalised_HU
            if (normalised_HU - normalised_flows[int(element)-1])*100/normalised_flows[int(element)-1] < -10:
                fid.write('Element underperfused {} in {} {} for cluster {}\n'.format(
                    element, branch_mapping[cluster]['Lobe'],branch_mapping[cluster]['branch_name'], cluster))
                fid.write('with baseline flow fraction of {} and signal fraction of {}\n'.format(
                    normalised_flows[int(element) - 1], normalised_HU))
                fid.write('signal to flow is less by {} percent.\n\n'.format(100*(normalised_HU - normalised_flows[int(element)-1])/normalised_flows[int(element)-1]))
