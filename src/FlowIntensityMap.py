import os
import cv2
import pydicom
import numpy as np

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

def writeExNodeFile(path, filename, coords):
    """
    Write out an ex Node file with the given coords.
    :param path: the path you want the file to be written at
    :param filename: Filename to write to.
    :param coords: List of coordinate lists.
    :return: None
    """
    with open(os.path.join(path, filename), 'w') as f:
        f.write(' Group name : MAC\n')
        f.write(' #Fields=%1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 0\n')
        f.write('   y.  Value index= 2, #Derivatives= 0\n')
        f.write('   z.  Value index= 3, #Derivatives= 0\n')
        for i in range(len(coords)):
            f.write(' Node:         %.4d\n' % (i + 1001))
            f.write('    %s\n' % (coords[i][0]))
            f.write('    %s\n' % (coords[i][1]))
            f.write('    %s\n' % (coords[i][2]))

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

    # # Try increasing voxel size if all intensities within voxel are zero (usually at edges)
    # multiplier = 0
    # while np.min(intensity) == 0 and multiplier < 10:
    #     multiplier += 1
    #     Voxel_xy = np.floor(spacing / x_y_res)
    #     Voxel_z = np.floor(spacing / z_res) * multiplier
    #     for i in range(len(coords)):
    #         if intensity[i] == 0:
    #             if coords[i, 2] < Voxel_z:
    #                 if coords[i, 2] <= 0:
    #                     nz = 0
    #                 elif coords[i, 2] > 0:
    #                     nz = coords[i, 2]
    #                 nnz = np.count_nonzero(Lung3D_double[(coords[i, 1] - Voxel_xy):(coords[i, 1] + Voxel_xy),
    #                                        (coords[i, 0] - Voxel_xy):(coords[i, 0] + Voxel_xy),
    #                                        : (Voxel_z + nz)])
    #                 if nnz == 0:
    #                     intensity[i] = 0
    #                 else:
    #                     intensity[i] = np.sum(Lung3D_double[(coords[i, 1] - Voxel_xy):(coords[i, 1] + Voxel_xy),
    #                                           (coords[i, 0] - Voxel_xy):(coords[i, 0] + Voxel_xy),
    #                                           : (Voxel_z + nz)]) / nnz
    #             elif coords[i, 2] > (len(Lung3D) - Voxel_z):
    #                 if coords[i, 2] >= len(Lung3D):
    #                     nz = 0
    #                 elif coords[i, 2] < len(Lung3D):
    #                     nz = len(Lung3D) - coords[i, 2]
    #                 nnz = np.count_nonzero(Lung3D_double[(coords[i, 1] - Voxel_xy):(coords[i, 1] + Voxel_xy),
    #                                        (coords[i, 0] - Voxel_xy):(coords[i, 0] + Voxel_xy),
    #                                        (len(Lung3D) - (Voxel_z + nz)):])
    #                 if nnz == 0:
    #                     intensity[i] = 0
    #                 else:
    #                     intensity[i] = np.sum(Lung3D_double[(coords[i, 1] - Voxel_xy):(coords[i, 1] + Voxel_xy),
    #                                           (coords[i, 0] - Voxel_xy):(coords[i, 0] + Voxel_xy),
    #                                           (len(Lung3D) - (Voxel_z + nz)):]) / nnz
    #             else:
    #                 nnz = np.count_nonzero(Lung3D_double[(coords[i, 1] - Voxel_xy):(coords[i, 1] + Voxel_xy),
    #                                        (coords[i, 0] - Voxel_xy):(coords[i, 0] + Voxel_xy),
    #                                        (coords[i, 2] - Voxel_z):(coords[i, 2] + Voxel_z)])
    #                 if nnz == 0:
    #                     intensity[i] = 0
    #                 else:
    #                     intensity[i] = np.sum(Lung3D_double[(coords[i, 1] - Voxel_xy):(coords[i, 1] + Voxel_xy),
    #                                           (coords[i, 0] - Voxel_xy):(coords[i, 0] + Voxel_xy),
    #                                           (coords[i, 2] - Voxel_z):(coords[i, 2] + Voxel_z)]) / nnz
    #
    # # Set 0 values to closest intensity (if increasing voxel size doesn't work)
    # counter = 0
    # for i in range(len(intensity)):
    #     if intensity[i] == 0:
    #         counter += 1
    #         intensity_sub1 = 0
    #         intensity_sub2 = 0
    #         k = i
    #         l = i
    #         while intensity_sub1 == 0 and intensity_sub2 == 0:
    #             k += 1
    #             l -= 1
    #             if k < len(intensity):
    #                 intensity_sub1 = intensity[k]
    #             elif l > 0:
    #                 intensity_sub2 = intensity[l]
    #         if intensity_sub1 != 0:
    #             intensity[i] = intensity_sub1
    #         else:
    #             intensity[i] = intensity_sub2

    intensity = intensity - 1024  # HU_offset
    return intensity

    # Output file
    # export_file = open(f'/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/{subject}_{lobe}_avg0g.exdata', 'w')
    # export_file.write(' Group name: colourful_spots\n')
    # export_file.write('#Fields=2\n')
    # export_file.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    # export_file.write('  x.  Value index= 1, #Derivatives=0\n')
    # export_file.write('  y.  Value index= 2, #Derivatives=0\n')
    # export_file.write('  z.  Value index= 3, #Derivatives=0\n')
    # export_file.write(' 2) UnknownField, field, rectangular cartesian, #Components=1\n')
    # export_file.write('  1.  Value index= 4, #Derivatives=0\n')
    # for j in range(len(intensity)):
    #     export_file.write(f' Node:      {j}\n')
    #     export_file.write(f'   {A.data[j, 2]}  {A.data[j, 3]}  {A.data[j, 4]}   {intensity[j]}\n')
    # export_file.close()


subject_name = 'Alfred1'
Lung_Path = os.path.join('/hpc/bsha219/lung/Data/CTEPH/', subject_name, 'Pre/scans/5-DE_CTPA_Fines__F_0.7/resources/DICOM')
Mask_Path = os.path.join('/hpc/bsha219/lung/Data/CTEPH/', subject_name, 'Pre/Lung/PTKLungMask')
lungs, images, masks = createDICOMStackAndMask(Lung_Path, Mask_Path, 0)
coords = read_coordinates('/hpc/bsha219/lung/Data/CTEPH/CTEPH10/FRC/Lobe/RLL_DataGrid.exdata')
intensity = voxel_avg_intensity(0.646484, 0.7, 5.2, coords, lungs)


