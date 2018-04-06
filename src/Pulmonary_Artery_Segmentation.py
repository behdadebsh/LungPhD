import sys
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('TkAgg')
# matplotlib.use('Qt4Agg')
# import matplotlib
# matplotlib.use('WXAgg')
from PIL import ImageTk, Image
import Tkinter as tk
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
# from mayavi import mlab
import numpy as np
from load_scan import load_scan, get_pixels_hu
from Lvl_Set_Alg import lvlset


# def click_roi(event, x, y, flags, params):
#     global roi_number, pts
#     while roi_number <= 3:
#         if event == cv2.EVENT_LBUTTONDOWN:
#             pts.append((x, y))
#             roi_number = roi_number + 1
#         elif event == cv2.EVENT_LBUTTONDBLCLK:
#             break
#     print pts
#     print pts[0, 0], pts[0, 1]
#     return pts
coords = []


# def callback(event):
#     print "clicked at", event.x, event.y


def onclick(event):
    # ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d'%(
        event.x, event.y)
    global coords
    coords.append((event.x, event.y))
    return coords

# def on_mouse(event, x, y, flags, params):
#     if event == cv2.EVENT_LBUTTONDOWN:
#         print 'Seed: ' + str(x) + ', ' + str(y)
#         global ix, iy
#         ix, iy = x, y


def main():
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    # data_dir = '/hpc/bsha219/lung/Data/ST12/Raw/DICOMS'
    patient_scans = load_scan('/hpc/bsha219/lung/Data/ST12/Raw/DICOMS')
    patient_images = get_pixels_hu(patient_scans)
    global artery_mask
    artery_mask = []
    for k in range(120, 170):
        global coords
        coords = []
        img = patient_images[k]
        for i in range(img.shape[0]):
            for j in range(img.shape[0]):
                if img[i][j] > 0:
                    img[i][j] *= 20
        blurred_f = ndimage.gaussian_filter(img, sigma=3)
        filter_blurred_f = ndimage.gaussian_filter(blurred_f, sigma=1)
        filter_blurred_f2 = filters.median_filter(filter_blurred_f, 3)
        alpha = 100
        sharpened = blurred_f + alpha * (blurred_f - filter_blurred_f2)

        root = tk.Tk()
        im = Image.fromarray(sharpened)
        photo = ImageTk.PhotoImage(im)
        canvas = tk.Canvas(root, width=512, height=512)
        canvas.pack()
        canvas.create_image(256, 256, image=photo)
        canvas.bind("<Button-1>", onclick)
        root.mainloop()

        # fig = plt.figure()
        # plt.gray()
        # plt.imshow(img, cmap='gray')
        # cid = fig.canvas.mpl_connect('button_press_event', onclick)
        # plt.show()
        # fig.canvas.mpl_disconnect(cid)
        # time.sleep(10)
        seed = np.zeros(img.shape)
        for jj in range(len(coords)):
            seed[int(coords[jj][1]), int(coords[jj][0])] = 1
        seg, phi, its = lvlset(sharpened, seed, max_its=1000, display=True, alpha=0.2)
        segment = seg * 1.0
        matplotlib.image.imsave('/hpc/bsha219/Python/Behdad/Out_arterial_mask/' + str(i) + '.png', segment)
        # fig.canvas.mpl_disconnect(cid)
        # plt.close()
        # artery_mask.append(segment)
    # artery_mask = np.asarray(artery_mask)
    # src = mlab.pipeline.scalar_field(artery_mask)
    # mlab.pipeline.volume(src, vmin=0, vmax=0.8)


if __name__ == '__main__':
    main()
