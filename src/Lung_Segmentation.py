#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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

@author: Behdad Shaarbaf Ebrahimi & Kaggle
            This code has been modified by Behdad Shaarbaf Ebrahimi (UoA - ABI) for academic purposes.
"""

import sys

import numpy as np
import dicom
import os
import scipy
import scipy.ndimage as ndimage
#from mayavi import mlab
import matplotlib.pyplot as plt
from skimage import measure, morphology, segmentation



# Load the scans in given folder path
def load_scan(path_files):
    slices = [dicom.read_file(path_files + '/' + s) for s in os.listdir(path_files)]
    slices.sort(key=lambda x: int(x.InstanceNumber))
    try:
        slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
    except:
        slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)
        
    for s in slices:
        s.SliceThickness = slice_thickness
        
    return slices


def get_pixels_hu(scans):
    image = np.stack([s.pixel_array for s in scans])
    # Convert to int16 (from sometimes int16), 
    # should be possible as values should always be low enough (<32k)
    image = image.astype(np.int16)

    # Set outside-of-scan pixels to 0
    # The intercept is usually -1024, so air is approximately 0
    image[image == -2000] = 0
    
    # Convert to Hounsfield units (HU)
    intercept = scans[0].RescaleIntercept
    slope = scans[0].RescaleSlope
    
    if slope != 1:
        image = slope * image.astype(np.float64)
        image = image.astype(np.int16)
        
    image += np.int16(intercept)
    
    return np.array(image, dtype=np.int16)


#sl = input('Which Slice do you want to mask? Enter it here: ')
#print ("Original Slice")
#plt.imshow(patient_images[sl], cmap='gray')
#plt.show()


# Some of the starting Code is taken from ArnavJain, since it's more readable then my own
def generate_markers(image):
    # Creation of the internal Marker
    marker_internal = image < -500
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
    left_lung = marker_internal_labels == 1
    right_lung = marker_internal_labels == 2
    marker_internal = marker_internal_labels > 0
    
    # Creation of the external Marker
    
    external_a = ndimage.binary_dilation(marker_internal, iterations=30)
    external_b = ndimage.binary_dilation(marker_internal, iterations=60)
    
    ex_left_a = ndimage.binary_dilation(left_lung, iterations=30)
    ex_left_b = ndimage.binary_dilation(left_lung, iterations=60)
    
    ex_right_a = ndimage.binary_dilation(right_lung, iterations=30)
    ex_right_b = ndimage.binary_dilation(right_lung, iterations=60)
    
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

# Show some example markers from the middle

# test_patient_internal, test_patient_external, test_patient_watershed, left_lung, right_lung, ex_left, ex_right, watershed_left, watershed_right = generate_markers(patient_images[slice])
# print ("Internal Marker")
# plt.imshow(test_patient_internal, cmap='gray')
# plt.show()
# print ("Left Lung Marker")
# plt.imshow(left_lung, cmap='gray')
# plt.show()
# print ("Right Lung Marker")
# plt.imshow(right_lung, cmap='gray')
# plt.show()


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
    # seg = np.where(seg != 0, image, -2000*np.ones((512, 512)))

    return seg, watershed, marker_internal, marker_external, marker_watershed, watershed_left, watershed_right


def resample(image, scan, new_spacing=[1, 1, 1]):
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


def vesselness(image):
    segmented, watershed, marker_int, marker_ext, marker_watershed, watershed_left, watershed_right = lung_segment(image)
    
    vessels = segmented > 0
    vessel_mask = ndimage.binary_opening(vessels, iterations=2)
    vessels = image * vessel_mask
#    print ("Vessel mask")
#    plt.imshow(vessel_mask, cmap='gray')
#    plt.show()
#    print ("Vessel Segmentation")
#    plt.imshow(vessels, cmap='gray')
#    plt.show()
    return vessel_mask, vessels
    

###############################################################################


def main():
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    patient_scans = load_scan(data_dir)
    patient_images = get_pixels_hu(patient_scans)
    imgs, spacing = resample(patient_images, patient_scans, [1, 1, 1])
    segmented = []
    wat = []
    vessel = []
    for i in range(len(imgs)):
        seg, watershed, mark_int, mark_ext, mark_watershed, wat_left, wat_right = lung_segment(imgs[i])
        ves_m, ves = vesselness(imgs[i])
        segmented.append(seg)
        wat.append(watershed)
        vessel.append(ves)
    segmented = np.asarray(segmented)
    wat = np.asarray(wat)
    vessel = np.asarray(vessel)
    wat[wat==128] = 0
    # src = mlab.pipeline.scalar_field(wat)
    #mlab.pipeline.volume(src, vmin=0, vmax=0.9)

    #mlab.pipeline.image_plane_widget(src,
                                #plane_orientation='x_axes',
                                #slice_index=10,
                                #)
    #mlab.pipeline.image_plane_widget(src,
                                #plane_orientation='y_axes',
                                #slice_index=10,
                                #)

if __name__ == '__main__':
    main()