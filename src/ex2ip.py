import sys
import re

file_name = sys.argv[1]

######################## READING EXNODE DATA ##############################

try:
    with open(f"{file_name}.exnode", "r") as exnode:
        file_data = exnode.readlines()
except FileNotFoundError:
    sys.exit(f"NO {file_name}.exnode EXISTS")

total_lines = len(file_data)
total_nodes = 0
old_node = []
index = {}
new_node = []
node = {}

line = 0
while line < total_lines:
    current_line = file_data[line].strip()
    if current_line.startswith("Node:"):
        total_nodes += 1
        words = current_line.split()
        if len(words) < 2:
            print(f"Warning: Line {line + 1} does not contain enough data: {current_line}")
            line += 1
            continue
        node_id = words[1]
        old_node.append(node_id)
        index[node_id] = total_nodes
        new_node.append(total_nodes)

        node[total_nodes] = []
        if line + 1 < total_lines:
            node[total_nodes].append(file_data[line + 1].strip())
        if line + 2 < total_lines:
            node[total_nodes].append(file_data[line + 2].strip())
        if line + 3 < total_lines:
            node[total_nodes].append(file_data[line + 3].strip())

        # print(f"Processed Node {total_nodes} (ID {node_id}): {node[total_nodes]}")

        line += 4  # Skip to the next node
    else:
        line += 1


#### Print statements for debugging
# print("Total nodes processed:", total_nodes)
# print("Nodes dictionary keys:", list(node.keys()))

######################## CREATING IPNODE FILE ##############################

try:
    with open(f"{file_name}.ipnode", "w") as ipnode:
        ipnode.write(" CMISS Version 2.0  ipnode File Version 2\n")
        ipnode.write(f" Heading: \n")
        ipnode.write(" \n")
        ipnode.write(f" The number of nodes is [   {total_nodes}]:    {total_nodes}\n")
        ipnode.write(" Number of coordinates [3]: 3\n")
        ipnode.write(" Do you want prompting for different versions of nj=1 [N]? Y\n")
        ipnode.write(" Do you want prompting for different versions of nj=2 [N]? Y\n")
        ipnode.write(" Do you want prompting for different versions of nj=3 [N]? Y\n")
        ipnode.write(" The number of derivatives for coordinate 1 is [0]: 0\n")
        ipnode.write(" The number of derivatives for coordinate 2 is [0]: 0\n")
        ipnode.write(" The number of derivatives for coordinate 3 is [0]: 0\n")
        ipnode.write(" \n")

        for n in range(1, total_nodes + 1):
            if n not in node:
                print(f"Warning: Node {n} not found in the data.")
                continue
            ipnode.write(f" Node number [    {n}]:     {n}\n")
            ipnode.write(" The number of versions for nj=1 is [1]:  1\n")
            ipnode.write(f" The Xj(1) coordinate is [ 0.00000E+00]:  {node[n][0] if len(node[n]) > 0 else ''}\n")
            ipnode.write(" The number of versions for nj=2 is [1]:  1\n")
            ipnode.write(f" The Xj(2) coordinate is [ 0.00000E+00]:  {node[n][1] if len(node[n]) > 1 else ''}\n")
            ipnode.write(" The number of versions for nj=3 is [1]:  1\n")
            ipnode.write(f" The Xj(3) coordinate is [ 0.00000E+00]:  {node[n][2] if len(node[n]) > 2 else ''}\n")
            ipnode.write(" \n")
except FileNotFoundError:
    sys.exit(f"NO {file_name}.ipnode EXISTS")

exelemfile = f"{file_name}.exelem"
ipelemfile = f"{file_name}.ipelem"

### Open exelem file ###
try:
    with open(exelemfile, "r") as exelem:
        exelem_lines = exelem.readlines()
except FileNotFoundError:
    sys.exit(f"\033[31mError: Can't open exelem file {exelemfile}\033[0m ")

### Reading "Group name"
line = exelem_lines[0]
if " Group name:" in line:
    group = re.search(r"\s*Group name:\s*(\w*)\s*", line).group(1)
else:
    sys.exit("\033[31mError: \"Group name\" not found in exelem header\033[0m ")

# Corrected line to count elements
NumberOfElems = sum(1 for line in exelem_lines if re.match(r"\s*Element:\s*\S*\s*0\s*0", line))
print(f"Total number of elements {NumberOfElems}")

### Exporting to ipelem ###
### Open ipelem file
try:
    with open(ipelemfile, "w") as ipelem, open(exelemfile, "r") as exelem:
        ### Write ipelem header
        ipelem.write(" CMISS Version 1.21 ipelem File Version 2\n")
        ipelem.write(f" Heading: {group}\n\n")
        ipelem.write(f" The number of elements is [{NumberOfElems:5d}]: {NumberOfElems:5d}\n")
        nj = 3

        nb = 1  # basis function number

        ### Loop over elems
        COUNT = 0

        for line in exelem:
            if re.match(r"\s*Element:\s*(\S*)\s*0\s*0", line):
                COUNT = int(re.search(r"\s*Element:\s*(\S*)\s*0\s*0", line).group(1))
                ipelem.write("\n")
                ipelem.write(f" Element number [{COUNT:5d}]: {COUNT:5d}\n")
                ipelem.write(f" The number of geometric Xj-coordinates is [{nj:1d}]: {nj:1d}\n")
                for i in range(1, nj + 1):
                    ipelem.write(f" The basis function type for geometric variable {i} is [{nb:1d}]: {nb:1d}\n")

                next(exelem)  # Skip the next line
                nodes = next(exelem).strip()
                ipelem.write(f" Enter the 2 global numbers for basis {nb}:       {nodes}\n")

except FileNotFoundError:
    sys.exit("\033[31mError: Can't open ipelem file\033[0m ")
