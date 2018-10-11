"""
Main Pulmonary Artery (MPA) segmentation example:

Make sure you have the required packages installed on you python interpreter/virtual environment before running this
script. Use the documentation for required packages on GitHub link: https://github.com/behdadebsh/LungPhD
Note: onclick function had to be repeated because the result of the click is saved in a global variable defined in
this function.

Author: Behdad Shaarbaf Ebrahimi
"""

from load_scan import load_scan, get_pixels_hu
from Lvl_Set_Alg import lvlset
import scipy
from PIL import ImageTk, Image
import Tkinter as tk
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import numpy as np


def onclick(event):
    # ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d'%(
        event.x, event.y)
    global coords
    coords.append((event.x, event.y))
    return coords


ct_images = load_scan('../examples/sample_inputs')
images_hu = get_pixels_hu(ct_images)
# global artery_mask
# artery_mask = []
for k in range(len(images_hu)):
    global coords
    coords = []
    img = images_hu[k]
    root = tk.Tk()
    im = Image.fromarray(img)
    photo = ImageTk.PhotoImage(im)
    canvas = tk.Canvas(root, width=512, height=512)
    canvas.pack()
    canvas.create_image(256, 256, image=photo)
    canvas.bind("<Button-1>", onclick)
    root.mainloop()
    for i in range(img.shape[0]):
        for j in range(img.shape[0]):
            if img[i][j] > 0:
                img[i][j] *= 12
    blurred_f = ndimage.gaussian_filter(img, sigma=3)
    filter_blurred_f = filters.median_filter(img, 3)
    filter_blurred_f2 = ndimage.gaussian_filter(filter_blurred_f, sigma=1)

    alpha = 10
    sharpened = filter_blurred_f + alpha * (filter_blurred_f - filter_blurred_f2)
    seed = np.zeros(img.shape)
    for jj in range(len(coords)):
        seed[int(coords[jj][1]), int(coords[jj][0])] = 1
    seg, phi, its = lvlset(sharpened, seed, max_its=1800, display=True, alpha=0.2, thresh=0.1)
    segment = seg * 1.0
    scipy.misc.imsave('./results/MPA%d.jpg' % (k+1), segment)
