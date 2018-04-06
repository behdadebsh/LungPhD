"""
Created on Fri Feb  2 15:54:59 2018

This script is for Lung Segmentation from DICOM images. Loading the scans from a directory and
transforms the pixels to Hounsfield Units. The features implemented in these codes are written
in a way to read the series of DICOM images from a folder and convert the voxel values into
hounsfield unit numbers. In order to segment the lung, two functions are made to perform this action.
First, markers are made via thresholding and then labeled. Labels are sorted from smaller areas to
bigger areas. Therefore, the 2 largest areas would be lungs. Second, using some morphological tools
and using the output from the first function to make lung mask and segment slices. Another feature of this
script is to visualise the outputs in 2D and 3D graphics.

@author: Behdad Shaarbaf Ebrahimi
    The functions are originated from Kaggle and modified by Behdad Shaarbaf Ebrahimi (UoA - ABI) for academic purposes.
"""

import sys, os
import numpy as np
from load_scan import load_scan, get_pixels_hu
import scipy
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('wx')
from mayavi import mlab
from skimage import measure, morphology, segmentation


def generate_markers(image):
    # Creation of the internal Marker
    marker_internal = image < -500  # Lung Tissue threshold
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
    left_lung = marker_internal_labels == 1  # marking left lung with 1 values
    right_lung = marker_internal_labels == 2  # marking right lung with 2 values
    marker_internal = marker_internal_labels > 0
    
    # Creation of the external Marker
    
    external_a = ndimage.binary_dilation(marker_internal, iterations=20)
    external_b = ndimage.binary_dilation(marker_internal, iterations=60)
    
    ex_left_a = ndimage.binary_dilation(left_lung, iterations=20)
    ex_left_b = ndimage.binary_dilation(left_lung, iterations=20)
    
    ex_right_a = ndimage.binary_dilation(right_lung, iterations=20)
    ex_right_b = ndimage.binary_dilation(right_lung, iterations=20)
    
    marker_external = external_b ^ external_a
    marker_ex_left = ex_left_b ^ ex_left_a
    marker_ex_right = ex_right_b ^ ex_right_a
    
    # Creation of the Watershed Marker matrix
    img_length = len(image[1])
    marker_watershed = np.zeros((img_length, img_length), dtype=np.int)
    watershed_left = np.zeros((img_length, img_length), dtype=np.int)
    watershed_right = np.zeros((img_length, img_length), dtype=np.int)
    
    marker_watershed += marker_internal * 255
    marker_watershed += marker_external * 128
    watershed_left += left_lung * 255
    watershed_left += marker_ex_left * 128
    watershed_right += right_lung * 255
    watershed_right += marker_ex_right * 128
    
    return marker_internal, marker_external, marker_watershed, left_lung, right_lung, marker_ex_left, marker_ex_right,\
        watershed_left, watershed_right


# Function using watershed algorithm ro to lung segmentation
def lung_segment(image):
    # Creation of the markers:
    marker_internal, marker_external, marker_watershed, left_lung, right_lung, ex_left, ex_right, watershed_left, watershed_right = generate_markers(image)
    
    # Creation of the Sobel-Gradient:
    sobel_filtered_dx = ndimage.sobel(image, 1)
    sobel_filtered_dy = ndimage.sobel(image, 0)
    sobel_gradient = np.hypot(sobel_filtered_dx, sobel_filtered_dy)
    sobel_gradient *= 255.0 / np.max(sobel_gradient)
    
    # Watershed algorithm:
    watershed = morphology.watershed(sobel_gradient, marker_watershed)
    watershed_left = morphology.watershed(sobel_gradient, watershed_left)
    watershed_right = morphology.watershed(sobel_gradient, watershed_right)
    img_length = len(image[1])
    for i in range(img_length):
        for j in range(img_length):
            if watershed_left[i, j] == 128:
                watershed_left[i, j] = 0
            if watershed_left[i, j] == 255:
                watershed_left[i, j] = 1
            if watershed_right[i, j] == 128:
                watershed_right[i, j] = 0
            if watershed_right[i, j] == 255:
                watershed_right[i, j] = 1
    
    seg_left = image * watershed_left
    seg_right = image * watershed_right
    seg = seg_left + seg_right
    # seg = np.where(seg != 0, image, -2000*np.ones((512, 512))) # optional line in order give every voxel except lung a value of -2000 hounsfield unit

    return seg, watershed, marker_internal, marker_external, marker_watershed, watershed_left, watershed_right


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
    patient_scans = load_scan('/hpc/bsha219/lung/Data/ST12/Raw/DICOMS')
    patient_images = get_pixels_hu(patient_scans)
    # imgs, spacing = downsample(patient_images, patient_scans, [1, 1, 1])
    segmented = []
    wat = []
    for i in range(len(patient_images)):
        seg, watershed, mark_int, mark_ext, mark_watershed, wat_left, wat_right = lung_segment(patient_images[i])
        watershed[watershed == 128] = 0
        watershed[watershed == 255] = 1
        # matplotlib.image.imsave('/hpc/bsha219/Python/Behdad/Lung_masks/' + str(i) + '.png', watershed)
        segmented.append(seg)
        wat.append(watershed)
    segmented = np.asarray(segmented)
    wat = np.asarray(wat)
    wat[wat == 128] = 0
    src1 = mlab.pipeline.scalar_field(wat)
    mlab.pipeline.volume(src1, vmin=0, vmax=0.8)
    # src2 = mlab.pipeline.scalar_field(segmented)
    # mlab.pipeline.volume(src2, vmin=0, vmax=0.8)
    # mlab.pipeline.image_plane_widget(src1,plane_orientation='x_axes',slice_index=10)
    # mlab.pipeline.image_plane_widget(src2,plane_orientation='y_axes',slice_index=10)


if __name__ == '__main__':
    main()
