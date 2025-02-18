#!/usr/bin/env python3
import sys, os
import numpy as np
from stl import mesh
import math

def main():
    # Load the STL file
    if len(sys.argv) > 1:
        mesh_file = str(sys.argv[1])

        if len(sys.argv) > 2:
            lung = str(sys.argv[2])
        else:
            lung = 'lung'
    else:
        root_folder = '/hpc/bsha219/lung/Data'
        study = 'CTEPH'
        protocol = 'Pre'
        subject = 'P-Z2R78'  # 'Alfred1'
        # data, header = nrrd.read(os.path.join(root_folder, study, subject, protocol, 'Vessel/MPA_Mask', subject + '_' + protocol + '.nrrd'))
        mesh_file = os.path.join(root_folder, study, subject, protocol, 'Lobe', 'Lobe mesh', 'LobeSurfaceMesh_RM.stl')
        #mesh_file = '/home/hcha184/Data/Tairawhiti_study/ants_test/004/004_lung.stl'
        # mesh_file = '/home/hcha184/Data/Alfred_CTEPH/Alfred05_translated/Post/Alfred05_Post_lung.stl'
        lung = 'lung'


    stl_file = mesh.Mesh.from_file(mesh_file)
    vertices = stl_file.vectors
    # flipped_vertices = vertices * np.array([-1, -1, 1]) # NEEDED to account for DICOM (LPS) NIFIT (RAS) convention
    flipped_vertices = vertices * np.array([1, 1, 1])  # No flip
    # Flatten the list of vertices
    vertices = flipped_vertices.reshape(-1, 3)

    num_vertices = len(vertices)
    if num_vertices > 2000000:
        skip = 100
    elif 1000000 < num_vertices >= 2000000:
        skip = 50
    else:
        skip = 15
    num_trimmed = math.floor(num_vertices / skip)

    output_exdata = mesh_file[:-3]+'exdata'
    if lung == 'Left':
        f = open(output_exdata, "w+")
        f.write("Group name: surface_Left\n")
        f.write("#Fields=1\n")
        f.write("1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
        f.write(" x.  Value index= 1, #Derivatives= 0\n")
        f.write(" y.  Value index= 2, #Derivatives= 0\n")
        f.write(" z.  Value index= 3, #Derivatives= 0\n")

        node_count = 10000
        for i in range(num_trimmed):
            node_count += 1
            num_temp = i*skip+1
            coord = [vertices[num_temp][0], vertices[num_temp][1], vertices[num_temp][2]]
            f.write("Node: %d\r\n" % node_count)
            f.write("   %3.5f  %3.5f  %3.5f\n" % (coord[0], coord[1], coord[2]))
        f.close()
    elif lung == 'Right':
        f = open(output_exdata, "w+")
        f.write("Group name: surface_Right\n")
        f.write("#Fields=1\n")
        f.write("1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
        f.write(" x.  Value index= 1, #Derivatives= 0\n")
        f.write(" y.  Value index= 2, #Derivatives= 0\n")
        f.write(" z.  Value index= 3, #Derivatives= 0\n")

        node_count = 200000
        for i in range(num_trimmed):
            node_count += 1
            num_temp = i * skip + 1
            coord = [vertices[num_temp][0], vertices[num_temp][1], vertices[num_temp][2]]
            f.write("Node: %d\r\n" % node_count)
            f.write("   %3.5f  %3.5f  %3.5f\n" % (coord[0], coord[1], coord[2]))
        f.close()

    else:
        f = open(output_exdata, "w+")
        f.write("Group name: surface_%s\n" % lung)
        f.write("#Fields=1\n")
        f.write("1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
        f.write(" x.  Value index= 1, #Derivatives= 0\n")
        f.write(" y.  Value index= 2, #Derivatives= 0\n")
        f.write(" z.  Value index= 3, #Derivatives= 0\n")

        node_count = 10000
        for i in range(num_trimmed):
            node_count += 1
            num_temp = i * skip + 1
            coord = [vertices[num_temp][0], vertices[num_temp][1], vertices[num_temp][2]]
            f.write("Node: %d\r\n" % node_count)
            f.write("   %3.5f  %3.5f  %3.5f\n" % (coord[0], coord[1], coord[2]))
        f.close()

if __name__ == "__main__":
    main()
