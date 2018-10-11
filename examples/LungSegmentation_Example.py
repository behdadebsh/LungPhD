"""
Lung segmentation running example:

Make sure you have the required packages installed on you python interpreter/virtual environment before running this
script. Use the documentation for required packages on GitHub link: https://github.com/behdadebsh/LungPhD


Author: Behdad Shaarbaf Ebrahimi
"""

from load_scan import load_scan, get_pixels_hu
from Lung_Segmentation import lung_segment
import matplotlib.pyplot as plt          # needed for viewing the images
from PIL import Image                    # needed for saving the images on hard drive
import numpy as np


ct_images = load_scan('../examples/sample_inputs')
images_hu = get_pixels_hu(ct_images)

for i in range(len(images_hu)):
    lungs, left_lung, right_lung = lung_segment(images_hu[i])
    lungs = np.uint8(lungs)
    binary_im = Image.fromarray(lungs)   # Using PIL here
    binary_im.save('./results/LungMask%d.jpg' % i+1, quality=100)  # Saving the slice into a file
    plt.figure()

    plt.subplot(2, 2, 1)
    plt.imshow(right_lung)

    plt.subplot(2, 2, 2)
    plt.imshow(left_lung)

    plt.subplot(2, 2, 3)
    plt.imshow(lungs)
    plt.show()
