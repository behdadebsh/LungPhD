from load_scan import load_scan, get_pixels_hu
from Lung_Segmentation import lung_segment, generate_markers
import scipy.ndimage as ndimage
import matplotlib
import matplotlib.pyplot as plt


patient_scans = load_scan('/hpc/bsha219/lung/Data/ST12/Raw/DICOMS')
patient_images = get_pixels_hu(patient_scans)
for i in range(len(patient_images)):
    marker_internal, marker_external, marker_watershed, left_lung, right_lung, ex_left, ex_right, watershed_left, watershed_right = generate_markers(patient_images[i])
    seg, watershed, mark_int, mark_ext, mark_watershed, wat_left, wat_right = lung_segment(patient_images[i])
    watershed[watershed == 128] = 0
    watershed[watershed == 255] = 1
    ex = watershed - marker_internal*1.0
    vessel = ndimage.binary_opening(ex)
    # matplotlib.image.imsave('/hpc/bsha219/Python/Behdad/Inside_vessels/' + str(i) + '.png', vessel)
# plt.imshow(vessel, cmap='gray')
# plt.show()
