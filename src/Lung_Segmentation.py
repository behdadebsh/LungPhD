"""

This script is for Lung Segmentation from DICOM images. Loading the scans from a directory and
transforms the pixels to Hounsfield Units. The features implemented in these codes are written
in a way to read the series of DICOM images from a folder and convert the voxel values into
hounsfield unit numbers. In order to segment the lung, two functions are made to perform this action.
First, markers are made via thresholding and then labeled. Labels are sorted from smaller areas to
bigger areas. Therefore, the 2 largest areas would be lungs. Second, using watershed algorithm
and using the output from the first function to make lung mask and segment slices. Another feature of this
script is to visualise the outputs in 2D and 3D graphics using mayavi.

@author: Behdad Shaarbaf Ebrahimi (UoA - ABI) for academic purposes.
"""

import sys
import numpy as np
from load_scan import load_scan, get_pixels_hu
import scipy
import scipy.ndimage as ndimage
from PIL import Image
# import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('wx')          # depending on your working backend you should use
# matplotlib.use('Qt4Agg')      # depending on your working backend you should use
# from mayavi import mlab
from skimage import measure, morphology, segmentation


def generate_markers(image):
    # Creation of the internal Marker
    marker_internal = image < -400  # Lung Tissue threshold
    marker_internal = segmentation.clear_border(marker_internal)
    marker_internal_labels = measure.label(marker_internal)
    areas = [r.area for r in measure.regionprops(marker_internal_labels)]
    areas.sort()
    if len(areas) > 2:
        for region in measure.regionprops(marker_internal_labels):
            if region.area < areas[-2]:
                for coordinates in region.coords:
                    marker_internal_labels[coordinates[0], coordinates[1]] = 0

    marker_internal_labels = measure.label(marker_internal_labels)

    # Creation of the external Marker
    
    external_a = ndimage.binary_dilation(marker_internal_labels, iterations=20)
    external_b = ndimage.binary_dilation(marker_internal_labels, iterations=60)
    marker_external = external_b ^ external_a

    # Creation of the Watershed Marker matrix
    img_length = len(image[1])
    marker_watershed = np.zeros((img_length, img_length), dtype=np.float64)

    marker_watershed += marker_internal_labels * 255
    marker_watershed += marker_external * 128

    return marker_internal_labels, marker_watershed


# Function using watershed algorithm ro to lung segmentation
def lung_segment(image):
    # Creation of the markers:
    marker_internal, marker_watershed = generate_markers(image)
    
    # Creation of the Sobel-Gradient:
    sobel_filtered_dx = ndimage.sobel(image, 1)
    sobel_filtered_dy = ndimage.sobel(image, 0)
    sobel_gradient = np.hypot(sobel_filtered_dx, sobel_filtered_dy)

    # Watershed algorithm:
    watershed = morphology.watershed(sobel_gradient, marker_watershed)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            if watershed[i, j] == 128:
                watershed[i, j] = 0
    lung_centroids = [l.centroid for l in measure.regionprops(watershed)]
    labels = [label.label for label in measure.regionprops(watershed)]
    if lung_centroids[0][1] > 0.5 * len(image[1]) and labels[0] == 255:
        left_lung = watershed == 255  # marking left lung
        right_lung = watershed == 510  # marking right lung
    else:
        left_lung = watershed == 510
        right_lung = watershed == 255
    left_lung[left_lung != 0] = 1
    right_lung[right_lung != 0] = 2
    left_lung = left_lung * 5.0
    right_lung = right_lung * 17.0
    lungs = right_lung + left_lung

    return lungs, left_lung, right_lung


# a function to downsample the image stack to make the code run faster
def downsample(image, scan, new_spacing=[1, 1, 1]):
    # Determine current pixel spacing
    spacing = map(float, ([scan[0].SliceThickness] + scan[0].PixelSpacing))
    spacing = np.array(list(spacing))

    resize_factor = spacing / new_spacing
    new_real_shape = image.shape * resize_factor
    new_shape = np.round(new_real_shape)
    real_resize_factor = new_shape / image.shape
    new_spacing = spacing / real_resize_factor

    image = scipy.ndimage.interpolation.zoom(image, real_resize_factor)

    return image, new_spacing


def main():
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    # data_dir = '/hpc/bsha219/lung/Data/ST12/Raw/DICOMS'
    # data_dir = '/hpc/bsha219/lung/Data/P2BRP257-H12076/FRC/Raw/DICOMS'
    patient_scans = load_scan('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Raw/DICOMS')
    patient_images = get_pixels_hu(patient_scans)
# imgs, spacing = downsample(patient_images, patient_scans, [1, 1, 1])
#     segmented = []

    for i in range(len(patient_images)):
        lungs, left_lung, right_lung = lung_segment(patient_images[i])
    # segmented.append(lungs)
    # segmented = np.asarray(segmented)
        lungs = np.uint8(lungs)
    # to save them as a stack of masks in a directory
        binary_im = Image.fromarray(lungs)
        binary_im.save('/hpc/bsha219/Python/Behdad/Lung_masks/LungMask%.4d.jpg' % i, quality=100)
    # scipy.misc.imsave('/hpc/bsha219/Python/Behdad/Lung_masks/LungMask%.4d.png' % i, lungs)

    # src1 = mlab.pipeline.scalar_field(segmented)
    # mlab.pipeline.volume(src1, vmin=0, vmax=0.8)
    # mlab.pipeline.image_plane_widget(src1,plane_orientation='x_axes',slice_index=10)


if __name__ == '__main__':
    main()
