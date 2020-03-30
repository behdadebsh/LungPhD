from load_scan import load_scan, get_pixels_hu
from Lung_Segmentation import lung_segment, generate_markers
import scipy.ndimage as ndimage
import matplotlib
# import matplotlib.pyplot as plt


patient_scans = load_scan('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Raw/DICOMS')
patient_images = get_pixels_hu(patient_scans)
for slice in range(len(patient_images)):
    marker_internal, marker_watershed = generate_markers(patient_images[slice])
    lungs, l_lung, r_lung = lung_segment(patient_images[slice])
    marker_watershed[marker_watershed == 128] = 0
    vessel = lungs - marker_internal
    for i in range(len(vessel[0])):
        for j in range(len(vessel[1])):
            if vessel[i][j] == 5 or vessel[i][j] == 17:
                vessel[i][j] = 1
            else:
                vessel[i][j] = 0
    vessel = ndimage.binary_opening(vessel)
    # plt.imshow(vessel)
    # plt.show()
    matplotlib.image.imsave('/hpc/bsha219/Python/Behdad/Inside_vessels/' + str(slice) + '.png', vessel)
