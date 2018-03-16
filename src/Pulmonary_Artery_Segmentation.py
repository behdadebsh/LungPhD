import sys
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import numpy as np
from load_scan import load_scan, get_pixels_hu
from Lvl_Set_Alg import lvlset

# def vesselness(image):
#     segmented, watershed, marker_int, marker_ext, marker_watershed, watershed_left, watershed_right = lung_segment(image)
#
#     vessels = segmented > 0
#     vessel_mask = ndimage.binary_opening(vessels, iterations=1)
#     vessels = image * vessel_mask
#     print ("Vessel mask")
#     plt.imshow(vessel_mask, cmap='gray')
#     plt.show()
#     print ("Vessel Segmentation")
#     plt.imshow(vessels, cmap='gray')
#     plt.show()
#     return vessel_mask, vessels

def main():
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    # data_dir = '/hpc/bsha219/lung/Data/ST12/Raw/DICOMS'
    patient_scans = load_scan('/hpc/bsha219/lung/Data/ST12/Raw/DICOMS')
    patient_images = get_pixels_hu(patient_scans)
    img = patient_images[150]
    for i in range(512):
        for j in range(512):
            if img[i][j] > 0:
                img[i][j] = img[i][j] * 20
    blurred_f = ndimage.gaussian_filter(img, sigma=3)
    filter_blurred_f = ndimage.gaussian_filter(blurred_f, sigma=1)
    alpha = 70
    sharpened = blurred_f + alpha * (blurred_f - filter_blurred_f)

    mask = np.zeros(img.shape)
    mask[300, 300] = 1
    # plt.imshow(img, cmap='gray')
    # plt.show()
    lvlset(sharpened, mask, max_its=2000, display=True, alpha=0.2)


if __name__ == '__main__':
    main()

