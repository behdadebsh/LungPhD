import sys
# import cv2
# import pyautogui
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import numpy as np
from load_scan import load_scan, get_pixels_hu
from Lvl_Set_Alg import lvlset
# import argparse


# initialize the list of reference points and boolean indicating
# whether cropping is being performed or not
# refPt = []
# cropping = False
#
#
# def click_and_crop(event, x, y, flags, param):
#     # grab references to the global variables
#     global refPt, cropping
#
#     # if the left mouse button was clicked, record the starting
#     # (x, y) coordinates and indicate that cropping is being
#     # performed
#     if event == cv2.EVENT_LBUTTONDOWN:
#         refPt = [(x, y)]
#         # cropping = True
#
#     # check to see if the left mouse button was released
#     # elif event == cv2.EVENT_LBUTTONUP:
#         # record the ending (x, y) coordinates and indicate that
#         # the cropping operation is finished
#         # refPt.append((x, y))
#         # cropping = False
#
#         # draw a rectangle around the region of interest
#         # cv2.rectangle(image, refPt[0], refPt[1], (0, 255, 0), 2)
#         # cv2.imshow("image", image)
#     return refPt


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
    img = img * 30
    # load the image, clone it, and setup the mouse callback function
    # image = img
    # clone = image.copy()
    # cv2.namedWindow("image")
    # cv2.setMouseCallback("image", click_and_crop)
    #
    # # keep looping until the 'q' key is pressed
    # while True:
    #     # display the image and wait for a keypress
    #     cv2.imshow("image", image)
    #     key = cv2.waitKey(1) & 0xFF
    #
    #     # if the 'r' key is pressed, reset the cropping region
    #     if key == ord("r"):
    #         image = clone.copy()
    #
    #     # if the 'c' key is pressed, break from the loop
    #     elif key == ord("c"):
    #         break
    #
    # # if there are two reference points, then crop the region of interest
    # # from teh image and display it
    # if len(refPt) == 2:
    #     roi = clone[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]]
    #     cv2.imshow("ROI", roi)
    #     cv2.waitKey(0)
    mask = np.zeros(img.shape)
    mask[300, 280] = 1
    # plt.imshow(img, cmap='gray')
    # plt.show()
    lvlset(img, mask, max_its=1000, display=True, alpha=1.0)


if __name__ == '__main__':
    main()

