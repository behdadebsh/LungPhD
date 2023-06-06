def data_node_lobes(Lobe_Path, intensity_lll, intensity_lul, intensity_rll, intensity_rml, intensity_rul, blood, air):
    # DATA POINT INTENSITIES
    # the file CTEPH4_RLL_avg.exdata (or for any lobe) contains the
    # seed point coords and the normalised intensities at these points.

    import os

    # RLL
    # os.chdir(Intensity_Path)
    # with open(Data_File_RLL, 'r') as fid_right:
    #     N = Data_Points_RLL  # # data points
    #     intensity_rll = [[0] * 2 for _ in range(N)]  # preallocate list
    #     intensity_rll[0][0] = int(fid_right.readline().split()[-1])
    #     intensity_rll[0][1] = float(fid_right.readline().split()[-1])
    #     for k in range(1, N):
    #         intensity_rll[k][0] = int(fid_right.readline().split()[-1])
    #         intensity_rll[k][1] = float(fid_right.readline().split()[-1])
    intensities_rll = intensity_rll

    # RML
    # os.chdir(Intensity_Path)
    # with open(Data_File_RML, 'r') as fid_right:
    #     N = Data_Points_RML  # # data points
    #     intensity_rml = [[0] * 2 for _ in range(N)]  # preallocate list
    #     intensity_rml[0][0] = int(fid_right.readline().split()[-1])
    #     intensity_rml[0][1] = float(fid_right.readline().split()[-1])
    #     for k in range(1, N):
    #         intensity_rml[k][0] = int(fid_right.readline().split()[-1])
    #         intensity_rml[k][1] = float(fid_right.readline().split()[-1])
    intensities_rml = intensity_rml

    # RUL
    # os.chdir(Intensity_Path)
    # with open(Data_File_RUL, 'r') as fid_right:
    #     N = Data_Points_RUL  # # data points
    #     intensity_rul = [[0] * 2 for _ in range(N)]  # preallocate list
    #     intensity_rul[0][0] = int(fid_right.readline().split()[-1])
    #     intensity_rul[0][1] = float(fid_right.readline().split()[-1])
    #     for k in range(1, N):
    #         intensity_rul[k][0] = int(fid_right.readline().split()[-1])
    #         intensity_rul[k][1] = float(fid_right.readline().split()[-1])
    intensities_rul = intensity_rul

    # LLL
    # os.chdir(Intensity_Path)
    # with open(Data_File_LLL, 'r') as fid_left:
    #     N = Data_Points_LLL  # # data points
    #     intensity_lll = [[0] * 2 for _ in range(N)]  # preallocate list
    #     intensity_lll[0][0] = int(fid_left.readline().split()[-1])
    #     intensity_lll[0][1] = float(fid_left.readline().split()[-1])
    #     for k in range(1, N):
    #         intensity_lll[k][0] = int(fid_left.readline().split()[-1])
    #         intensity_lll[k][1] = float(fid_left.readline().split()[-1])
    intensities_lll = intensity_lll

    # LUL
    # os.chdir(Intensity_Path)
    # with open(Data_File_LUL, 'r') as fid_left:
    #     N = Data_Points_LUL  # # data points
    #     intensity_lul = [[0] * 2 for _ in range(N)]  # preallocate list
    #     intensity_lul[0][0] = int(fid_left.readline().split()[-1])
    #     intensity_lul[0][1] = float(fid_left.readline().split()[-1])
    #     for k in range(1, N):
    #         intensity_lul[k][0] = int(fid_left.readline().split()[-1])
    #         intensity_lul[k][1] = float(fid_left.readline().split()[-1])
    intensities_lul = intensity_lul

    import numpy as np

    # Define file paths
    os.chdir(Lobe_Path)
    # lobe_path = "path/to/lobe/folder"
    filename1 = "RLL_mapping.txt"
    filename2 = "RML_mapping.txt"
    filename3 = "RUL_mapping.txt"
    filename4 = "LLL_mapping.txt"
    filename5 = "LUL_mapping.txt"

    # Import mapping data
    # Read and process file 1
    with open(f"{Lobe_Path}/{filename1}", 'r') as file:
        lines = file.readlines()

    A = np.zeros((len(lines), 2), dtype=int)
    for i, line in enumerate(lines):
        values = line.split()
        A[i, 0] = int(values[0])
        A[i, 1] = int(values[1])

    # Read and process file 2
    with open(f"{Lobe_Path}/{filename2}", 'r') as file:
        lines = file.readlines()

    B = np.zeros((len(lines), 2), dtype=int)
    for i, line in enumerate(lines):
        values = line.split()
        B[i, 0] = int(values[0])
        B[i, 1] = int(values[1])

    # Read and process file 3
    with open(f"{Lobe_Path}/{filename3}", 'r') as file:
        lines = file.readlines()

    C = np.zeros((len(lines), 2), dtype=int)
    for i, line in enumerate(lines):
        values = line.split()
        C[i, 0] = int(values[0])
        C[i, 1] = int(values[1])

    # Read and process file 4
    with open(f"{Lobe_Path}/{filename4}", 'r') as file:
        lines = file.readlines()

    D = np.zeros((len(lines), 2), dtype=int)
    for i, line in enumerate(lines):
        values = line.split()
        D[i, 0] = int(values[0])
        D[i, 1] = int(values[1])

    # Read and process file 5
    with open(f"{Lobe_Path}/{filename5}", 'r') as file:
        lines = file.readlines()

    E = np.zeros((len(lines), 2), dtype=int)
    for i, line in enumerate(lines):
        values = line.split()
        E[i, 0] = int(values[0])
        E[i, 1] = int(values[1])

    # Concatenate all matrices
    # nodal_intensities = np.concatenate((A, B, C, D, E), axis=0)

    l = len(A) + len(B) + len(C) + len(D) + len(E)
    nodal_intensities = np.zeros((l, 2))
    # # Map nodal intensities to corresponding lung regions
    # intensities_rll = np.loadtxt("intensities_rll.txt")
    # intensities_rml = np.loadtxt("intensities_rml.txt")
    # intensities_rul = np.loadtxt("intensities_rul.txt")
    # intensities_lll = np.loadtxt("intensities_lll.txt")
    # intensities_lul = np.loadtxt("intensities_lul.txt")

    for k in range(len(A)):
        nodal_intensities[k, 0] = int(A[k, 1])
        nodal_intensities[k, 1] = intensities_rll[k]

    for k in range(len(A), len(A) + len(B)):
        k2 = k - len(A)
        nodal_intensities[k, 0] = int(B[k2, 1])
        nodal_intensities[k, 1] = intensities_rml[k2]

    for k in range(len(A) + len(B), len(A) + len(B) + len(C)):
        k2 = k - (len(A) + len(B))
        nodal_intensities[k, 0] = int(C[k2, 1])
        nodal_intensities[k, 1] = intensities_rul[k2]

    for k in range(len(A) + len(B) + len(C), len(A) + len(B) + len(C) + len(D)):
        k2 = k - (len(A) + len(B) + len(C))
        nodal_intensities[k, 0] = int(D[k2, 1])
        nodal_intensities[k, 1] = intensities_lll[k2]

    for k in range(len(A) + len(B) + len(C) + len(D), l):
        k2 = k - (len(A) + len(B) + len(C) + len(D))
        nodal_intensities[k, 0] = int(E[k2, 1])
        nodal_intensities[k, 1] = intensities_lul[k2]

    # Normalize intensities
    range_ = blood - air

    for k in range(l):
        if nodal_intensities[k, 1] < air:
            nodal_intensities[k, 1] = air
        elif nodal_intensities[k, 1] > blood:
            nodal_intensities[k, 1] = blood
        nodal_intensities[k, 1] = (nodal_intensities[k, 1] - air) / range_
    return nodal_intensities

