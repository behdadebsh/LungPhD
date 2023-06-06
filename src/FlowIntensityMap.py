import os
import cv2
import pydicom
import numpy as np
from global_nodes import global_nodes
from data_node_lobes import data_node_lobes

def createDICOMStackAndMask(Lung_Path, Mask_Path, option):
    """
    :param Lung_Path: Path to the DICOM images directory
    :param Mask_Path: Path to the mask images directory
    :param option: Selecting the lung. 0 stands for both, 1 stands for left, and 2 for right lung
    :return: A list of masked image arrays. Each list component is a 2D array.
    """
    if option not in [1, 2]:
        option = 0  # For both lungs selection

    # Selection for which lung to import and mask.
    # (0)Both lungs (1)Left lung (2)Right lung
    Lung_Selection = option

    # Intensity values of each lung in the mask. Output from PASS makes
    # left lung = 5 and right lung = 17
    Right_Lung_Intensity = 17
    Left_Lung_Intensity = 5

    # Read DICOM files and sort them based on InstanceNumber
    Lung_ImagesDir = sorted(os.listdir(Lung_Path), key=lambda f: int(pydicom.dcmread(os.path.join(Lung_Path, f)).InstanceNumber))
    # Mask_ImagesDir = sorted([f for f in os.listdir(Mask_Path) if f.endswith('.tiff')], key=lambda f: int(f.split('.')[0]))
    Mask_ImagesDir = sorted([f for f in os.listdir(Mask_Path) if f.endswith('.tiff')], key=lambda f: f.split('.')[0])
    # Mask_ImagesDir = [f for f in os.listdir(Mask_Path) if f.endswith('.tiff')]

    NumOfImages = len(Lung_ImagesDir)
    Lung3D = []  # Create list to hold all lung images
    Mask3D = []  # Create list to hold all mask images
    Image3D = []  # Create list to hold all the images

    # Reverse the order of mask images
    # Mask_ImagesDir = Mask_ImagesDir[::-1]

    # Loop through lung images dir and store segmented images
    for k in range(2, NumOfImages):
        # Import current Lung and Mask images to work on
        Lung_CurrentImage = pydicom.dcmread(os.path.join(Lung_Path, Lung_ImagesDir[k])).pixel_array
        Mask_CurrentImage = cv2.imread(os.path.join(Mask_Path, Mask_ImagesDir[k - 2]), cv2.IMREAD_GRAYSCALE)

        if Lung_Selection == 0:
            # Both lungs.
            # Threshold mask image to binary.
            _, Mask_CurrentImage = cv2.threshold(Mask_CurrentImage, 0, 255, cv2.THRESH_BINARY)
        elif Lung_Selection == 1:
            # Left lung
            Mask_CurrentImage = cv2.inRange(Mask_CurrentImage, Left_Lung_Intensity, Left_Lung_Intensity)
        else:
            # Right lung
            Mask_CurrentImage = cv2.inRange(Mask_CurrentImage, Right_Lung_Intensity, Right_Lung_Intensity)

        # Flip the mask image vertically
        Mask_CurrentImage = np.flip(Mask_CurrentImage, axis=0)

        # Remove noise from the edges of the lungs using erosion filter
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
        Mask_CurrentImage = cv2.erode(Mask_CurrentImage, kernel)

        # Multiply the mask and the lung image to acquire the Region of Interest
        Masked_Image = cv2.bitwise_and(Lung_CurrentImage, Lung_CurrentImage, mask=Mask_CurrentImage)

        # Store the masked image in the list
        Lung3D.append((Masked_Image))
        Image3D.append(Lung_CurrentImage)
        Mask3D.append(Mask_CurrentImage)

    return Lung3D, Image3D, Mask3D

def read_coordinates(filename):
    coordinates = []
    with open(filename) as f_in:
        lines = (line.rstrip() for line in f_in)
        file_content = list(line for line in lines if line)
    for file_line in range(len(file_content)):
        if (file_content[file_line].split()[0]) == 'Node:':
            x = float(file_content[file_line + 1].split()[0])
            y = float(file_content[file_line + 1].split()[1])
            z = float(file_content[file_line + 1].split()[2])
            coordinates.append((x, y, z))
    return np.array(coordinates)

def writeExNodeFileWithField(path, filename, coords, field):
    """
    Write out an ex Node file with the given coords.
    :param path: the path you want the file to be written at
    :param filename: Filename to write to.
    :param coords: List of coordinate lists.
    :return: None
    """
    with open(os.path.join(path, filename + '.exnode'), 'w') as f:
        f.write(' Group name : voxel_intensity\n')
        f.write(' #Fields=2\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 0\n')
        f.write('   y.  Value index= 2, #Derivatives= 0\n')
        f.write('   z.  Value index= 3, #Derivatives= 0\n')
        f.write(' 2) avg_intensity, field, rectangular cartesian, #Components=1\n')
        f.write('  1.  Value index= 4, #Derivatives=0\n')
        # for j in range(len(intensity)):
        #     f.write(f' Node:      {j}\n')
        #     f.write(f'   {coords[j, 2]}  {coords[j, 3]}  {coords[j, 4]}   {intensity[j]}\n')
        # export_file.close()
        for i in range(len(coords)):
            f.write(' Node:         %.4d\n' % (i + 1001))
            f.write('    %s\n' % (coords[i][0]))
            f.write('    %s\n' % (coords[i][1]))
            f.write('    %s\n' % (coords[i][2]))
            f.write('    %s\n' % (field[i]))
        f.close()

def voxel_avg_intensity(x_y_res,z_res,spacing,coords,Lung3D):

    #x_y_res =  # Set the value of x_y_res
    #z_res =  # Set the value of z_res
    #spacing =  # Set the value of spacing

    coords[:, 0:2] = np.round(
        (coords[:, 0:2] - x_y_res / 2) / x_y_res)  # coordinate transformation from cartesian to pixel
    coords[:, 2] = -1 * np.round(
        (coords[:, 2] - z_res / 2) / z_res)  # coordinate transformation from cartesian to image slice

    # Find intensities at specified locations, averaging algorithm here*******
    Voxel_xy = np.floor((spacing / x_y_res) / 2)  # In theory half the number of pixels between each seed point in every direction (rounded down)
    Voxel_z = np.floor((spacing / z_res) / 2)
    pixels = Lung3D[0].shape
    Lung3D_double = np.zeros((pixels[0], pixels[1], len(Lung3D)))
    for i in range(len(Lung3D)):
        Lung3D_double[:, :, i] = Lung3D[i][:, :]
    Lung3D_double = np.array(Lung3D_double)

    intensity = np.zeros(len(coords))

    for i in range(len(coords)):
        if coords[i, 2] < Voxel_z:
            if coords[i, 2] <= 0:
                nz = 0
            elif coords[i, 2] > 0:
                nz = coords[i, 2]
            nnz = np.count_nonzero(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                   int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                   : int(Voxel_z + nz)])
            if nnz == 0:
                intensity[i] = 0
            else:
                intensity[i] = np.sum(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                      int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                      : int(Voxel_z + nz)]) / nnz
        elif coords[i, 2] > (len(Lung3D) - Voxel_z):
            if coords[i, 2] >= len(Lung3D):
                nz = 0
            elif coords[i, 2] < len(Lung3D):
                nz = len(Lung3D) - coords[i, 2]
            nnz = np.count_nonzero(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                   int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                   int(len(Lung3D) - (Voxel_z + nz)):])
            if nnz == 0:
                intensity[i] = 0
            else:
                intensity[i] = np.sum(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                      int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                      int(len(Lung3D) - (Voxel_z + nz)):]) / nnz
        else:
            nnz = np.count_nonzero(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                   int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                   int(coords[i, 2] - Voxel_z):int(coords[i, 2] + Voxel_z)])
            if nnz == 0:
                intensity[i] = 0
            else:
                intensity[i] = np.sum(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                      int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                      int(coords[i, 2] - Voxel_z):int(coords[i, 2] + Voxel_z)]) / nnz

    # Try increasing voxel size if all intensities within voxel are zero (usually at edges)
    multiplier = 0
    while np.min(intensity) == 0 and multiplier < 10:
        multiplier += 1
        Voxel_xy = np.floor(spacing / x_y_res)
        Voxel_z = np.floor(spacing / z_res) * multiplier
        for i in range(len(coords)):
            if intensity[i] == 0:
                if int(coords[i, 2]) < Voxel_z:
                    if int(coords[i, 2]) <= 0:
                        nz = 0
                    elif int(coords[i, 2]) > 0:
                        nz = int(coords[i, 2])
                    nnz = np.count_nonzero(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                           int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                           : (Voxel_z + nz)])
                    if nnz == 0:
                        intensity[i] = 0
                    else:
                        intensity[i] = np.sum(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                              int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                              : (Voxel_z + nz)]) / nnz
                elif int(coords[i, 2]) > (len(Lung3D) - Voxel_z):
                    if int(coords[i, 2]) >= len(Lung3D):
                        nz = 0
                    elif int(coords[i, 2]) < len(Lung3D):
                        nz = len(Lung3D) - coords[i, 2]
                    nnz = np.count_nonzero(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                           int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                           (len(Lung3D) - (Voxel_z + nz)):])
                    if nnz == 0:
                        intensity[i] = 0
                    else:
                        intensity[i] = np.sum(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                              int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                              (len(Lung3D) - (Voxel_z + nz)):]) / nnz
                else:
                    nnz = np.count_nonzero(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                           int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                           int(coords[i, 2] - Voxel_z):int(coords[i, 2] + Voxel_z)])
                    if nnz == 0:
                        intensity[i] = 0
                    else:
                        intensity[i] = np.sum(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                              int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                              int(coords[i, 2] - Voxel_z):int(coords[i, 2] + Voxel_z)]) / nnz

    # Set 0 values to closest intensity (if increasing voxel size doesn't work)
    counter = 0
    for i in range(len(intensity)):
        if intensity[i] == 0:
            counter += 1
            intensity_sub1 = 0
            intensity_sub2 = 0
            k = i
            l = i
            while intensity_sub1 == 0 and intensity_sub2 == 0:
                k += 1
                l -= 1
                if k < len(intensity):
                    intensity_sub1 = intensity[k]
                elif l > 0:
                    intensity_sub2 = intensity[l]
            if intensity_sub1 != 0:
                intensity[i] = intensity_sub1
            else:
                intensity[i] = intensity_sub2

    intensity = intensity - 1024  # HU_offset
    return intensity

def get_num_elements(artery_file):
    with open(artery_file) as f_in:
        lines = (line.rstrip() for line in f_in)
        file_content = list(line for line in lines if line)
    for file_line in (file_content):
        if (file_line.split()[0]) == 'Element:':
            num_elems = int(file_line.split()[1])
        # else:
        #     print("Artery file corrupt.")
        #     exit(0)
    return num_elems

def intensity_mapping_avg(subject_path, subject_name, num_elems, intensity_rll, intensity_rml, intensity_rul, intensity_lll, intensity_lul, blood, air):
    # vessel_path = os.path.join(subject_path, 'Vessel')
    # artery_file = vessel_path + '/' + subject_name + '_Artery_Full.exnode'
    # num_elems = get_num_elements(artery_file)
    elem_global_nodes, elem_gen, elem_children = global_nodes(vessel_path, artery_file, num_elems)
    Lobe_Path = os.path.join(subject_path, 'Lobe')
    Intensity_Path = os.path.join(subject_path, 'Intensity_mapping')
    isExist = os.path.exists(Intensity_Path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(Intensity_Path)
    # Data_Points_LLL = len(intensity_lll)
    # Data_Points_LUL = len(intensity_lul)
    # Data_Points_RLL = len(intensity_rll)
    # Data_Points_RML = len(intensity_rml)
    # Data_Points_RUL = len(intensity_rul)
    nodal_intensities = data_node_lobes(Lobe_Path, intensity_lll, intensity_lul, intensity_rll,
                                        intensity_rml, intensity_rul, blood, air)

    # Calculate the sum of intensities
    Intensity_sum = np.sum(nodal_intensities[:, 1])

    # Calculate nodal intensity fractions
    nodal_intensity_fractions = nodal_intensities[:, 1] / Intensity_sum

    start = np.max(elem_gen)
    l = len(elem_gen)
    elem_intensities = np.zeros((l, 1))  # Adjust the size of elem_intensities

    for n in range(len(nodal_intensities)):
        elem_intensities[int(nodal_intensities[n, 0])-1] = nodal_intensity_fractions[n]

    for k in range(start, 0, -1):
        for m in range(l):
            generation = elem_gen[m]
            if generation == k:
                if elem_children[m, 0] == 0:
                    for n in range(len(nodal_intensities)):
                        if nodal_intensities[n, 0] == elem_global_nodes[m, 1]:
                            elem_intensities[m] = nodal_intensity_fractions[n]
                elif elem_children[m, 1] == 0:
                    elem_intensities[m] = elem_intensities[elem_children[m, 0]]
                else:
                    elem_intensities[m] = elem_intensities[elem_children[m, 0] - 1] + elem_intensities[elem_children[m, 1] - 1]

    # Export elem_intensities to a file
    export_file_1 = open(
        '/hpc/bsha219/lung/Data/CTEPH/' + subject_name + '/FRC/Intensity_mapping/' + subject_name + '_elem_intensity_map_python.exelem',
        'w')

    export_file_1.write(' Group name: flow_diff1\n')
    export_file_1.write(
        ' Shape.  Dimension=1\n #Scale factor sets=0\n #Nodes=0\n #Fields=1\n 1) flow_diff_fraction, field, rectangular cartesian, #Components=1\n   flow_diff_fraction.  l.Lagrange, no modify, grid based.\n   #xi1=1\n')

    for j in range(len(elem_intensities)):
        export_file_1.write(' Element:     ' + str(j+1) + ' ' + str(0) + ' ' + str(0) + '\n   Values:\n        ' + str(
            elem_intensities[j]).strip('[]') + '    ' + str(elem_intensities[j]).strip('[]') + '\n')

    export_file_1.close()


subject_name = 'Alfred1'
subject_path = '/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC'
Lung_Path = os.path.join('/hpc/bsha219/lung/Data/CTEPH/', subject_name, 'Pre/scans/5-DE_CTPA_Fines__F_0.7/resources/DICOM')
Mask_Path = os.path.join('/hpc/bsha219/lung/Data/CTEPH/', subject_name, 'Pre/Lung/PTKLungMask')
lungs, images, masks = createDICOMStackAndMask(Lung_Path, Mask_Path, 0)
coords_RLL = read_coordinates('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/RLL_DataGrid.exdata')
coords_RML = read_coordinates('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/RML_DataGrid.exdata')
coords_RUL = read_coordinates('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/RUL_DataGrid.exdata')
coords_LLL = read_coordinates('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/LLL_DataGrid.exdata')
coords_LUL = read_coordinates('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/LUL_DataGrid.exdata')
intensity_RLL = voxel_avg_intensity(0.646484, 0.7, 5.2, coords_RLL, lungs)
intensity_RML = voxel_avg_intensity(0.646484, 0.7, 5.2, coords_RML, lungs)
intensity_RUL = voxel_avg_intensity(0.646484, 0.7, 5.2, coords_RUL, lungs)
intensity_LLL = voxel_avg_intensity(0.646484, 0.7, 5.2, coords_LLL, lungs)
intensity_LUL = voxel_avg_intensity(0.646484, 0.7, 5.2, coords_LUL, lungs)
writeExNodeFileWithField('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe', 'test_RLL', coords_RLL, intensity_RLL)
writeExNodeFileWithField('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe', 'test_RUL', coords_RUL, intensity_RUL)
writeExNodeFileWithField('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe', 'test_RML', coords_RML, intensity_RML)
writeExNodeFileWithField('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe', 'test_LLL', coords_LLL, intensity_LLL)
writeExNodeFileWithField('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe', 'test_LUL', coords_LUL, intensity_LUL)
vessel_path = os.path.join(subject_path, 'Vessel')
artery_file = vessel_path + '/' + 'CTEPH10_Artery_Full.exelem'
num_elems = get_num_elements(artery_file)
intensity_mapping_avg(subject_path, 'CTEPH10', num_elems, intensity_RLL, intensity_RML, intensity_RUL, intensity_LLL, intensity_LUL, 417, -960)




