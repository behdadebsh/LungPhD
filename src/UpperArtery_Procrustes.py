
import numpy as np
from procrustes import procrustes
from utils import writeExNodeFile, writeipNodeFile

template_dir = '/hpc/bsha219/lung/GeometricModels/3D_Digitise/Templates/art_template.exnode'
centreline_dir = '/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/SurfaceFEMesh/centre_points.exdata'

with open(template_dir, 'rt') as f:
    template = f.read()
temp_lines = template.split('\n')
with open(centreline_dir, 'rt') as f:
    centreline = f.read()
centre_lines = centreline.split('\n')

# Reading Template Nodes of Trunk Vessel Coordinates

x = []  # x coordinates of template nodes
y = []  # y coordinates of template nodes
z = []  # z coordinates of template nodes
temp_coords = []
for line in range(len(temp_lines)):
    if 'Node:' in temp_lines[line]:
        x.append(float(temp_lines[line + 1]))
        y.append(float(temp_lines[line + 2]))
        z.append(float(temp_lines[line + 3]))
for j in range(11):
    temp_coords.append([x[j], y[j], z[j]])
temp_coords = np.asarray(temp_coords)

# Reading all the Template Nodes Coordinates

all_temp_coords = []
for j in range(len(x)):
    all_temp_coords.append([x[j], y[j], z[j]])
all_temp_coords = np.asarray(all_temp_coords)

# Reading Centreline Nodes Coordinates

xc = []  # x coordinates of centreline nodes
yc = []  # y coordinates of centreline nodes
zc = []  # z coordinates of centreline nodes
for line in range(len(centre_lines)):
    if 'Node:' in centre_lines[line]:
        line = line + 1
        xyz = centre_lines[line].split()
        xc.append(float(xyz[0]))
        yc.append(float(xyz[1]))
        zc.append(float(xyz[2]))


# Applying desired sequence for centreline points to correspond with template nodes

c_coords = []
node_seq = [1, 3, 5, 7, 27, 15, 23, 26, 11, 18, 21]
for row_number in node_seq:
    c_coords.append([xc[row_number - 1], yc[row_number - 1], zc[row_number - 1]])
c_coords = np.asarray(c_coords)


# Calculation of procrustes transformation on template to centreline

d, Z, tform = procrustes(c_coords, temp_coords)


# Applyting the rotation and translation on the Template coordinates

new_temp_coords = ((np.asmatrix(tform["rotation"]) * np.asmatrix(all_temp_coords).transpose()).transpose()) * tform["scale"] + tform["translation"]
a = new_temp_coords.tolist()


# Mapping Node #1 of template on Node #1 of Centreline
a[0][0] = xc[0]
a[0][1] = yc[0]
a[0][2] = zc[0]

# Mapping Node #2 of template on Node #3 of Centreline
a[1][0] = xc[2]
a[1][1] = yc[2]
a[1][2] = zc[2]

# Mapping Node #3 of template on Node #5 of Centreline
a[2][0] = xc[4]
a[2][1] = yc[4]
a[2][2] = zc[4]

# Mapping Node #4 of template on Node #7 of Centreline
a[3][0] = xc[6]
a[3][1] = yc[6]
a[3][2] = zc[6]

# Mapping Node #5 of template on Node #27 of Centreline
a[4][0] = xc[26]
a[4][1] = yc[26]
a[4][2] = zc[26]

# Mapping Node #6 of template on Node #15 of Centreline
a[5][0] = xc[14]
a[5][1] = yc[14]
a[5][2] = zc[14]

# Mapping Node #7 of template on Node #23 of Centreline
a[6][0] = xc[22]
a[6][1] = yc[22]
a[6][2] = zc[22]

# Mapping Node #8 and lower generation points relatively to Node #26
v821_x = xc[25] - a[7][0]
v821_y = yc[25] - a[7][1]
v821_z = zc[25] - a[7][2]
v821 = [v821_x, v821_y, v821_z]
right = [8, 12, 13, 14, 15, 18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33, 42, 43, 44, 45, 46, 47, 48, 49, 62, 63, 64, 65]
for i in right:
    a[i-1][0] += v821_x
    a[i-1][1] += v821_y
    a[i-1][2] += v821_z

# Mapping Node #9 of template on Node #11 of Centreline
a[8][0] = xc[10]
a[8][1] = yc[10]
a[8][2] = zc[10]

# Mapping Node #10 of template on Node #23 of Centreline
a[9][0] = xc[17]
a[9][1] = yc[17]
a[9][2] = zc[17]

# Mapping Node #11 and lower generation points relatively to Node #21
v1126_x = xc[20] - a[10][0]
v1126_y = yc[20] - a[10][1]
v1126_z = zc[20] - a[10][2]
v1126 = [v1126_x, v1126_y, v821_z]
left = [11, 16, 17, 22, 23, 24, 25, 34, 35, 36, 37, 38, 39, 40, 41, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69]
for i in left:
    a[i-1][0] += v1126_x
    a[i-1][1] += v1126_y
    a[i-1][2] += v1126_z

writeExNodeFile('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/SurfaceFEMesh/new_template.exnode', a)
writeipNodeFile('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/SurfaceFEMesh/new_template.ipnode', a)
