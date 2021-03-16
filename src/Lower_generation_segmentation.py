from load_scan import load_scan, get_pixels_hu
from Lung_Segmentation import generate_markers
import matplotlib.image as mpimg
import glob
import numpy as np
from PIL import Image
# import matplotlib.pyplot as plt


patient_scans = load_scan('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Raw/DICOMS')
patient_images = get_pixels_hu(patient_scans)
filelist_lung = glob.glob('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Lung/PTKLungMask/*.tif')
filelist_lung.sort()
lung_masks = np.array([np.array(Image.open(fname)) for fname in filelist_lung])
lung_masks[lung_masks != 0] = 1

for slice in range(len(patient_images)):
    marker_internal, marker_watershed = generate_markers(patient_images[slice])
    marker_watershed[marker_watershed == 255] = 1
    marker_watershed[marker_watershed == 510] = 1
    marker_watershed[marker_watershed != 1] = 0
    A = marker_watershed + lung_masks[slice]
    A = A*lung_masks[slice]
    A[A == 2] = 0
    mpimg.imsave('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/InsideVessels/Vessel%.4d.jpg' % slice, A)
