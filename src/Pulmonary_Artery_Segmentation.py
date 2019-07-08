import sys
import matplotlib
import matplotlib.pyplot as plt
import scipy
from PIL import ImageTk, Image
import Tkinter as tk
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import numpy as np
from load_scan import load_scan, get_pixels_hu
from Lvl_Set_Alg import lvlset


coords = []


def onclick(event):
    # ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d'%(
        event.x, event.y)
    global coords
    coords.append((event.x, event.y))
    return coords


def main():
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    # data_dir = '/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Raw/DICOMS'
    patient_scans = load_scan('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Raw/DICOMS')
    patient_images = get_pixels_hu(patient_scans)
    global artery_mask
    artery_mask = []
    for k in range(115, 215):
        global coords
        coords = []
        img = patient_images[k]
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
        # blurred_f = ndimage.gaussian_filter(img, sigma=3)
        filter_blurred_f = filters.median_filter(img, 3)
        filter_blurred_f2 = ndimage.gaussian_filter(filter_blurred_f, sigma=1)

        alpha = 12
        sharpened = filter_blurred_f + alpha * (filter_blurred_f - filter_blurred_f2)
        seed = np.zeros(img.shape)
        for jj in range(len(coords)):
            seed[int(coords[jj][1]), int(coords[jj][0])] = 1
        seg, phi, its = lvlset(sharpened, seed, max_its=1000, display=True, alpha=0.4)
        segment = seg * 1.0
        # scipy.misc.imsave('/hpc/bsha219/Python/Behdad/Out_arterial_mask/LungMask%.4d.jpg' % k, segment) # binary files
        # matplotlib.image.imsave('/hpc/bsha219/Python/Behdad/Out_arterial_mask/' + str(k) + '.png', segment)  # RGB files
        segment = np.uint8(segment)
        # to save them as a stack of masks in a directory
        binary_im = Image.fromarray(segment)
        binary_im.save('/hpc/bsha219/Python/Behdad/Lung_masks/LungMask%.4d.jpg' % k, quality=100)


if __name__ == '__main__':
    main()
