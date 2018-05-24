import numpy as np

#
# """
#
# This file reads surface data points and calculates centerlines and writes the coordinates as '*.exdata' format
# that is compatible with cmgui.
#
#
# """

file_dir = '/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/SurfaceFEMesh/surface_datapoints.exdata'
with open(file_dir, 'rt') as f:
    text = f.read()

start_lines = [7, 9, 167, 169, 39, 41, 103, 105, 71, 73, 135, 137]
lines = text.split('\n')
x = []
y = []
z = []
r = []

for line_nu in range(len(start_lines)):
    counter = 1
    tempx = []
    tempy = []
    tempz = []
    while counter < 9:
        st_line = lines[start_lines[line_nu]].split()
        tempx.append(float(st_line[0]))
        tempy.append(float(st_line[1]))
        tempz.append(float(st_line[2]))
        start_lines[line_nu] = start_lines[line_nu] + 4
        counter = counter + 1
        if counter == 8:
            x.append(float(np.mean(tempx)))
            y.append(float(np.mean(tempy)))
            z.append(float(np.mean(tempz)))
            r.append(np.sqrt((x[line_nu]-np.array(tempx[:]))**2 + (y[line_nu]-np.array(tempy[:]))**2 + (z[line_nu]-np.array(tempz[:]))**2))

radius = []
for j in range(len(r)):
    radius.append(np.mean(r[j]))


f1 = open("bif_centerline.exdata", 'w')
f1.write(' Group name : data points\n')
f1.write(' #Fields=2\n')
f1.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
f1.write('\tx.  Value index= 1, #Derivatives= 0\n')
f1.write('\ty.  Value index= 2, #Derivatives= 0\n')
f1.write('\tz.  Value index= 3, #Derivatives= 0\n')
f1.write(' 2) radius, field, rectangular cartesian, #Components=1\n')
f1.write('\tr.  Value index= 4, #Derivatives= 0\n')
f1.close()

with open('bif_centerline.exdata', 'at') as f:
    for i in range(1, 13):
        f.write(' Node:     %.4d\n' % i)
        f.write('   %s  %s  %s\n' % (x[i-1], y[i-1], z[i-1]))
        f.write('   %s\n' % radius[i-1])
