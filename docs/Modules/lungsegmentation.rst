====================
Lung_Segmentation.py
====================

This script is for Lung Segmentation from DICOM images. Loading the scans from a directory and
transforms the pixels to Hounsfield Units. The features implemented in these codes are written
in a way to read the series of DICOM images from a folder and convert the voxel values into
hounsfield unit numbers. In order to segment the lung, two functions are made to perform this action.
First, markers are made via thresholding and then labeled. Labels are sorted from smaller areas to
bigger areas. Therefore, the 2 largest areas would be lungs. Second, using watershed algorithm
and using the output from the first function to make lung mask and segment slices. Another feature of this script is to visualise the outputs in 2D and 3D graphics using mayavi. This script is adopted to separate one lung from another.

NOTES:

- The use of mayavi is optional. That is why it is commented out initially. In order to visualise the segmentation results in Python environment `mayavi <https://docs.enthought.com/mayavi/mayavi/>`_ package is exploited.

