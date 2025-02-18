import os
import cv2
import pydicom
import numpy as np
from global_nodes import global_nodes
from data_node_lobes import data_node_lobes
import pyvista as pv
from medpy.io import load
import SimpleITK as sitk
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cluster import KMeans
import time
from collections import defaultdict


def main():
    start = time.process_time()
    root_folder = '/hpc/bsha219/lung/Data'
    study = 'CTEPH'
    protocol = 'Pre'
    subject_name = 'Alfred1'
    masked_name = subject_name + '_' + protocol + '_masked.mha'
    #masked_name = subject_name + '_' + protocol + '_thresholdmasked_0.mha'
    print("For", subject_name, protocol)
    subject_path = os.path.join(root_folder, study, subject_name)
    masked_directory = os.path.join(subject_path, protocol, 'Lung')
    lobe_directory = os.path.join(subject_path, protocol, 'Lobe')
    ply_directory = os.path.join(subject_path, protocol, 'Lobe', 'Lobe_mesh')
    # Lung_Path = os.path.join(subject_path, protocol, 'Raw/DICOM')
    #Lung_Path = os.path.join(subject_path, 'Registration')
    # Mask_Path = os.path.join(subject_path, protocol, 'Raw/PTKLungMask')
    # lungs, images, masks = createDICOMStackAndMask(Lung_Path, Mask_Path, 0)

    # Read in masked metaimage
    img = sitk.ReadImage(os.path.join(masked_directory, masked_name))
    metadata = {}  # Metadata for masked image
    metadata['size'] = img.GetSize()
    metadata['spacing'] = img.GetSpacing()
    metadata['origin'] = img.GetOrigin()
    metadata['direction'] = img.GetDirection()

    # Convert the image to a numpy array
    lungs = sitk.GetArrayFromImage(img)

    if metadata['direction'][8] == 1: # Flip z for RAI to match RAS
        lungs = np.flip(lungs, axis=0) # ITK has array [z,y,x]
    # Reshape the array into a 3D array
    # img_array = img_array.reshape(metadata['size'][::-1])

    # Read lobe datapoints exdata and surface mesh ply
    coords_RLL = read_coordinates(os.path.join(lobe_directory, 'RLL_datapoints.exdata'))
    coords_RML = read_coordinates(os.path.join(lobe_directory, 'RML_datapoints.exdata'))
    coords_RUL = read_coordinates(os.path.join(lobe_directory, 'RUL_datapoints.exdata'))
    coords_LLL = read_coordinates(os.path.join(lobe_directory, 'LLL_datapoints.exdata'))
    coords_LUL = read_coordinates(os.path.join(lobe_directory, 'LUL_datapoints.exdata'))
    RUL_mesh = pv.read(os.path.join(ply_directory, 'RUL_surf.ply'))
    RML_mesh = pv.read(os.path.join(ply_directory, 'RML_surf.ply'))
    RLL_mesh = pv.read(os.path.join(ply_directory, 'RLL_surf.ply'))
    LUL_mesh = pv.read(os.path.join(ply_directory, 'LUL_surf.ply'))
    LLL_mesh = pv.read(os.path.join(ply_directory, 'LLL_surf.ply'))
    RUL_volume = RUL_mesh.volume
    RML_volume = RML_mesh.volume
    RLL_volume = RLL_mesh.volume
    LUL_volume = LUL_mesh.volume
    LLL_volume = LLL_mesh.volume
    RUL_spacing = (RUL_volume/6400)**(1/3)
    RML_spacing = (RML_volume/3200)**(1/3)
    RLL_spacing = (RLL_volume/8000)**(1/3)
    LUL_spacing = (LUL_volume/6400)**(1/3)
    LLL_spacing = (LLL_volume/8000)**(1/3)
    xy_res = metadata['spacing'][0]
    z_res = metadata['spacing'][2]

    direction = img.GetDirection()
    if len(direction) == 9:  # 3D case
        direction_matrix = np.reshape(direction, (3, 3))
    translation_vector = np.array(metadata['origin'])
    coords_RLL = coords_RLL - translation_vector
    coords_LLL = coords_LLL - translation_vector
    coords_RML = coords_RML - translation_vector
    coords_RUL = coords_RUL - translation_vector
    coords_LUL = coords_LUL - translation_vector

    # Map voxel intensity to datapoints
    print("=========== Started Mapping =============")
    coords_RLL, intensity_RLL = voxel_avg_intensity(xy_res, z_res, direction_matrix, RLL_spacing, coords_RLL, lungs)
    print("=========== RLL Mapped =============")
    coords_RML, intensity_RML = voxel_avg_intensity(xy_res, z_res, direction_matrix, RML_spacing, coords_RML, lungs)
    coords_RUL, intensity_RUL = voxel_avg_intensity(xy_res, z_res, direction_matrix, RUL_spacing, coords_RUL, lungs)
    coords_LLL, intensity_LLL = voxel_avg_intensity(xy_res, z_res, direction_matrix, LLL_spacing, coords_LLL, lungs)
    coords_LUL, intensity_LUL = voxel_avg_intensity(xy_res, z_res, direction_matrix, LUL_spacing, coords_LUL, lungs)

    coords_RLL = coords_RLL + translation_vector
    coords_LLL = coords_LLL + translation_vector
    coords_RML = coords_RML + translation_vector
    coords_RUL = coords_RUL + translation_vector
    coords_LUL = coords_LUL + translation_vector

    writeExNodeFileWithField(lobe_directory, 'RLL_datapoints', coords_RLL, intensity_RLL)
    writeExNodeFileWithField(lobe_directory, 'RUL_datapoints', coords_RUL, intensity_RUL)
    writeExNodeFileWithField(lobe_directory, 'RML_datapoints', coords_RML, intensity_RML)
    writeExNodeFileWithField(lobe_directory, 'LLL_datapoints', coords_LLL, intensity_LLL)
    writeExNodeFileWithField(lobe_directory, 'LUL_datapoints', coords_LUL, intensity_LUL)


    # Clustering based on HU values
    print("Started Clustering")
    # LLL = np.column_stack((coords_LLL, intensity_LLL))
    # LUL = np.column_stack((coords_LUL, intensity_LUL))
    # RLL = np.column_stack((coords_RLL, intensity_RLL))
    # RML = np.column_stack((coords_RML, intensity_RML))
    # RUL = np.column_stack((coords_RUL, intensity_RUL))

    scaler = MinMaxScaler()

    lobe = ['RUL', 'RML', 'RLL', 'LUL', 'LLL']
    for j in range(len(lobe)):
        print("Elem clustering for", subject_name, protocol, lobe[j])
        if lobe[j] == 'RUL':
            n_clusters = 6
            cluster_label = 1
            coords_norm = scaler.fit_transform(coords_RUL)
            orig_coords = coords_RUL
            intensity_norm = scaler.fit_transform(intensity_RUL.reshape(-1, 1))
            normalised_data = np.column_stack((coords_norm, intensity_norm))
        elif lobe[j] == 'RML':
            n_clusters = 4
            cluster_label = 7
            coords_norm = scaler.fit_transform(coords_RML)
            orig_coords = coords_RML
            intensity_norm = scaler.fit_transform(intensity_RML.reshape(-1, 1))
            normalised_data = np.column_stack((coords_norm, intensity_norm))
        elif lobe[j] == 'RLL':
            n_clusters = 6
            cluster_label = 11
            coords_norm = scaler.fit_transform(coords_RLL)
            orig_coords = coords_RLL
            intensity_norm = scaler.fit_transform(intensity_RLL.reshape(-1, 1))
            normalised_data = np.column_stack((coords_norm, intensity_norm))
        elif lobe[j] == 'LUL':
            n_clusters = 8
            cluster_label = 23
            coords_norm = scaler.fit_transform(coords_LUL)
            orig_coords = coords_LUL
            intensity_norm = scaler.fit_transform(intensity_LUL.reshape(-1, 1))
            normalised_data = np.column_stack((coords_norm, intensity_norm))
        else: # lobe = LLL
            n_clusters = 6
            cluster_label = 17
            coords_norm = scaler.fit_transform(coords_LLL)
            orig_coords = coords_LLL
            intensity_norm = scaler.fit_transform(intensity_LLL.reshape(-1, 1))
            normalised_data = np.column_stack((coords_norm, intensity_norm))

        kmeans = KMeans(n_clusters=n_clusters).fit(normalised_data)
        label_kmeans = kmeans.labels_
        clusters_norm = kmeans.cluster_centers_
        original_data = scaler.inverse_transform(normalised_data)
        original_centroids = scaler.inverse_transform(clusters_norm)
        clusters = original_centroids
        print(clusters)

        termdict = []

        for i in range(len(original_data)):
            termdict.append(dict(node_num=i + 1,
                                 x=orig_coords[i][0],
                                 y=orig_coords[i][1],
                                 z=orig_coords[i][2],
                                 intensity=original_data[i][3],
                                 cluster=label_kmeans[i] + cluster_label
                                 ))

        termdict = sorted(termdict, key=lambda k: k['node_num'])

        filename = lobe[j] + '_cluster_map_terminals_TTTTTT.exnode'
        export_path = os.path.join(subject_path, protocol, 'Intensity_mapping') + '/' + filename
        with open(export_path, 'w') as f:
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

            for i in range(len(original_data)):
                f.write(' Node:           %d\n' % (termdict[i]['node_num']))
                f.write('     %s\n' % float(termdict[i]['x']))
                f.write('     %s\n' % float(termdict[i]['y']))
                f.write('     %s\n' % float(termdict[i]['z']))
                f.write('     %s\n' % (termdict[i]['intensity']))
                f.write('     %s\n' % (termdict[i]['cluster']))

        clusters = defaultdict(list)
        for point in termdict:
            cluster_value = point['cluster']
            clusters[cluster_value].append([point['x'], point['y'], point['z']])
            clustered_points = dict(clusters)
        for cluster, points in clustered_points.items():
            writeipNodeFile(os.path.join(subject_path, protocol, 'Intensity_mapping', 'Cl_' + str(cluster) + '.ipdata'), points)
            writeExDataFile(os.path.join(subject_path, protocol, 'Intensity_mapping', 'Cl_' + str(cluster) + '.exdata'),
                            clusters, cluster)

    end = time.process_time()
    print("elapsed time:", end - start)


def writeExDataFile(filename, clusters, cluster_id, mean_radius_field=None):
    """
    Write out an ex data file for a specific cluster.
    :param filename: Filename to write to.
    :param clusters: defaultdict containing cluster data as {cluster_id: list of [x, y, z] coordinates}.
    :param cluster_id: The cluster ID to write out.
    :param mean_radius_field: Additional mean radius field to write out (optional).
    :return: None
    """
    coords = clusters[cluster_id]  # Extract coordinates for the given cluster ID

    with open(filename, 'w') as f:
        f.write(' Group name : cluster %d\n' % (cluster_id))
        f.write(' #Fields={0}\n'.format(1 if mean_radius_field is None else 2))
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('\tx.  Value index= 1, #Derivatives= 0\n')
        f.write('\ty.  Value index= 2, #Derivatives= 0\n')
        f.write('\tz.  Value index= 3, #Derivatives= 0\n')
        if mean_radius_field is not None:
            f.write(' 2) radius, field, rectangular cartesian, #Components=1\n')
            f.write('\tr.  Value index= 4, #Derivatives= 0\n')
        for i, coord in enumerate(coords):
            f.write(' Node:     %d\n' % (i + 1))
            f.write('   %s  %s  %s\n' % (coord[0], coord[1], coord[2]))
            if mean_radius_field is not None:
                f.write('   %s\n' % mean_radius_field[i])


def writeipNodeFile(filename, coords):
    """
    Write out an ipnode file with the given coords.
    :param filename: Filename to write to.
    :param coords: List of coordinate lists.
    :return: None
    """
    with open(filename, 'w') as f:
        f.write(' converted from exdata\n')
        for i in range(len(coords)):
            f.write(f"{i+1} {coords[i][0]:0.4E} {coords[i][1]:0.4E} {coords[i][2]:0.4E} 1.0 1.0 1.0\n")



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
    Mask_ImagesDir = sorted([f for f in os.listdir(Mask_Path) if f.endswith('.tiff')], key=lambda f: f.split('.')[0])
    image_data, image_header = load(os.path.join(Lung_Path))
    #image = sitk.ReadImage(os.path.join(Lung_Path))
    #image_array = sitk.GetArrayFromImage(image)
    # NumOfImages = len(image_array)
    NumOfImages = len(Lung_ImagesDir)
    Lung3D = []  # Create list to hold all lung images
    Mask3D = []  # Create list to hold all mask images
    Image3D = []  # Create list to hold all the images

    # Reverse the order of mask images
    # Mask_ImagesDir = Mask_ImagesDir[::-1]

    # Loop through lung images dir and store segmented images
    for k in range(NumOfImages):
        # Import current Lung and Mask images to work on
        #Lung_CurrentImage = image_data[k-2]
        Lung_CurrentImage = pydicom.dcmread(os.path.join(Lung_Path, Lung_ImagesDir[k])).pixel_array
        Mask_CurrentImage = cv2.imread(os.path.join(Mask_Path, Mask_ImagesDir[k]), cv2.IMREAD_GRAYSCALE)

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
        for i in range(len(coords)):
            f.write(' Node:         %.4d\n' % (i + 1001))
            f.write('    %s\n' % (coords[i][0]))
            f.write('    %s\n' % (coords[i][1]))
            f.write('    %s\n' % (coords[i][2]))
            f.write('    %s\n' % (field[i]))
        f.close()


def voxel_avg_intensity(x_y_res, z_res, direction_matrix, spacing, coords, Lung3D):
    # NOTE: code is this function assumes RAS image orientation
    coords = np.round((coords / np.array([x_y_res, x_y_res, z_res])).dot(np.linalg.inv(direction_matrix)))
    # coords[:,0:2] = np.round((coords[:,0:2]-x_y_res/2)/x_y_res)  # coordinate transformation from cartesian to pixel
    # coords[:,2] = -1 * np.round((coords[:,2]-z_res/2)/z_res)  # coordinate transformation from cartesian to image slice

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

            # test = Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
            #                        int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
            #                        int(coords[i, 2] - Voxel_z):int(coords[i, 2] + Voxel_z)]
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
                                           : int(Voxel_z + nz)])
                    if nnz == 0:
                        intensity[i] = 0
                    else:
                        intensity[i] = np.sum(Lung3D_double[int(coords[i, 1] - Voxel_xy):int(coords[i, 1] + Voxel_xy),
                                              int(coords[i, 0] - Voxel_xy):int(coords[i, 0] + Voxel_xy),
                                              : int(Voxel_z + nz)]) / nnz
                elif int(coords[i, 2]) > (len(Lung3D) - Voxel_z):
                    if int(coords[i, 2]) >= len(Lung3D):
                        nz = 0
                    elif int(coords[i, 2]) < len(Lung3D):
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

    # Set 0 values to closest intensity (if increasing voxel size doesn't work)
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

    #  TRANSFORM COORDS TO MIRROR AND MATCH ORIENTATION ###
    # transformation_matrix = np.array([[x_y_res, 0, 0],
    #                                   [0, x_y_res, 0],
    #                                   [0, 0, z_res]])
    coords = np.dot(coords*np.array([x_y_res, x_y_res, z_res]), direction_matrix)

    # intensity = intensity - 1024  # HU_offset
    return coords, intensity


def get_num_elements(artery_file):
    with open(artery_file) as f_in:
        lines = (line.rstrip() for line in f_in)
        file_content = list(line for line in lines if line)
    for file_line in (file_content):
        if (file_line.split()[0]) == 'Element:':
            num_elems = int(file_line.split()[1])
    return num_elems


def intensity_mapping_avg(subject_path, protocol, subject_name, vessel_path, artery_file, num_elems, intensity_rll, intensity_rml, intensity_rul, intensity_lll, intensity_lul, blood, air):

    elem_global_nodes, elem_gen, elem_children = global_nodes(vessel_path, artery_file, num_elems)
    Lobe_Path = os.path.join(subject_path, protocol, 'Lobe')
    Intensity_Path = os.path.join(subject_path, protocol, 'Intensity_mapping')
    isExist = os.path.exists(Intensity_Path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(Intensity_Path)
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
                    elem_intensities[m] = elem_intensities[elem_children[m, 0]-1]
                else:
                    elem_intensities[m] = elem_intensities[elem_children[m, 0]-1] + elem_intensities[elem_children[m, 1]-1]

    # Export elem_intensities to a file
    export_file_1 = open(Intensity_Path + '/' + subject_name + '_elem_intensity_map.exelem', 'w')

    export_file_1.write(' Group name: flow_diff1\n')
    export_file_1.write(
        ' Shape.  Dimension=1\n #Scale factor sets=0\n #Nodes=0\n #Fields=1\n 1) flow_diff_fraction, field, rectangular cartesian, #Components=1\n   flow_diff_fraction.  l.Lagrange, no modify, grid based.\n   #xi1=1\n')

    for j in range(len(elem_intensities)):
        export_file_1.write(' Element:     ' + str(j+1) + ' ' + str(0) + ' ' + str(0) + '\n   Values:\n        ' + str(
            elem_intensities[j]).strip('[]') + '    ' + str(elem_intensities[j]).strip('[]') + '\n')

    export_file_1.close()


def flow_comparison(subject_path, protocol, subject_name, num_elems):
    intensity_path = os.path.join(subject_path, protocol) + '/Intensity_mapping/'
    file = subject_name + '_elem_intensity_map.exelem'
    fid = open(intensity_path + file)
    intensity_fractions = np.zeros((num_elems, 1))  # preallocate array
    lines = (line.rstrip() for line in fid)
    int_file_content = list(line for line in lines if line)
    i = 0
    for file_line in range(len(int_file_content)):
        if (int_file_content[file_line].split()[0]) == 'Element:':
            intensity_fractions[i] = float(int_file_content[file_line + 2].split()[0])
            i = i + 1
    fid.close()

    # STEADY_FLOW_SOLUTION
    fid = open(os.path.join(subject_path, protocol, 'Perfusion', 'RM0.0') + '/' + subject_name + '_flow_perf.exelem')
    flow_fractions = np.zeros((num_elems, 1))  # preallocate array
    lines = (line.rstrip() for line in fid)
    flow_file_content = list(line for line in lines if line)
    i = 0
    for file_line in range(len(int_file_content)):
        if (flow_file_content[file_line].split()[0]) == 'Element:':
            flow_fractions[i] = float(flow_file_content[file_line + 2].split()[0])
            i = i + 1

    # Here comparison is ((B) - (A))/(B) i.e. difference as a fraction of (B)
    B_norm = flow_fractions[0]
    B = [elem / B_norm for elem in flow_fractions]  # normalize to input flow
    C = [(B[i] - intensity_fractions[i]) / B[i] for i in range(len(intensity_fractions))]  # difference as a fraction of flow solution

    # Export
    export_file_1 = open(intensity_path + '/' + subject_name + '_flow_diff_fractions.exelem',
        'w')
    export_file_1.write(' Group name: flow_diff1\n')
    export_file_1.write(
        ' Shape.  Dimension=1\n #Scale factor sets=0\n #Nodes=0\n #Fields=1\n 1) flow_diff_fraction, field, rectangular cartesian, #Components=1\n   flow_diff_fraction.  l.Lagrange, no modify, grid based.\n   #xi1=1\n')

    for j in range(len(C)):
        export_file_1.write(
            ' Element:     ' + str(j+1) + ' ' + str(0) + ' ' + str(0) + '\n   Values:\n        ' + str(
                C[j]).strip('[]') + '    ' + str(C[j]).strip('[]') + '\n')
    export_file_1.close()

if __name__ == '__main__':
    main()





