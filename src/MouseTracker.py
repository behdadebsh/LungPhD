import Tkinter as tk
from PIL import Image, ImageTk
import numpy as np
import matplotlib.pyplot as plt
import sys
from load_scan import load_scan, get_pixels_hu
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import glob
from LevelSet3D import chanvese3d
# import dateutil.parser


class MousePositionTracker(tk.Frame):
    """ Tkinter Canvas mouse position widget. """

    def __init__(self, canvas):
        self.canvas = canvas
        self.canv_width = self.canvas.cget('width')
        self.canv_height = self.canvas.cget('height')
        self.reset()

        # Create canvas cross-hair lines.
        xhair_opts = dict(dash=(3, 2), fill='white', state=tk.HIDDEN)
        self.lines = (self.canvas.create_line(0, 0, 0, self.canv_height, **xhair_opts),
                      self.canvas.create_line(0, 0, self.canv_width,  0, **xhair_opts))

    def cur_selection(self):
        return (self.start, self.end)

    def begin(self, event):
        self.hide()
        self.start = (event.x, event.y)  # Remember position (no drawing).

    def update(self, event):
        self.end = (event.x, event.y)
        self._update(event)
        self._command(self.start, (event.x, event.y))  # User callback.

    def _update(self, event):
        # Update cross-hair lines.
        self.canvas.coords(self.lines[0], event.x, 0, event.x, self.canv_height)
        self.canvas.coords(self.lines[1], 0, event.y, self.canv_width, event.y)
        self.show()

    def reset(self):
        self.start = self.end = None

    def hide(self):
        self.canvas.itemconfigure(self.lines[0], state=tk.HIDDEN)
        self.canvas.itemconfigure(self.lines[1], state=tk.HIDDEN)

    def show(self):
        self.canvas.itemconfigure(self.lines[0], state=tk.NORMAL)
        self.canvas.itemconfigure(self.lines[1], state=tk.NORMAL)

    def autodraw(self, command=lambda *args: None):
        """Setup automatic drawing; supports command option"""
        self.reset()
        self._command = command
        self.canvas.bind("<Button-1>", self.begin)
        self.canvas.bind("<B1-Motion>", self.update)
        self.canvas.bind("<ButtonRelease-1>", self.quit)

    def quit(self, event):
        self.hide()  # Hide cross-hairs.
        self.reset()


class SelectionObject:
    """ Widget to display a rectangular area on given canvas defined by two points
        representing its diagonal.
    """
    def __init__(self, canvas, select_opts):
        # Create attributes needed to display selection.
        self.canvas = canvas
        self.select_opts1 = select_opts
        self.width = self.canvas.cget('width')
        self.height = self.canvas.cget('height')

        # Options for areas outside rectanglar selection.
        select_opts1 = self.select_opts1.copy()  # Avoid modifying passed argument.
        select_opts1.update(state=tk.HIDDEN)  # Hide initially.
        # Separate options for area inside rectanglar selection.
        select_opts2 = dict(dash=(2, 2), fill='', outline='white', state=tk.HIDDEN)

        # Initial extrema of inner and outer rectangles.
        imin_x, imin_y,  imax_x, imax_y = 0, 0,  1, 1
        omin_x, omin_y,  omax_x, omax_y = 0, 0,  self.width, self.height

        self.rects = (
            # Area *outside* selection (inner) rectangle.
            self.canvas.create_rectangle(omin_x, omin_y,  omax_x, imin_y, **select_opts1),
            self.canvas.create_rectangle(omin_x, imin_y,  imin_x, imax_y, **select_opts1),
            self.canvas.create_rectangle(imax_x, imin_y,  omax_x, imax_y, **select_opts1),
            self.canvas.create_rectangle(omin_x, imax_y,  omax_x, omax_y, **select_opts1),
            # Inner rectangle.
            self.canvas.create_rectangle(imin_x, imin_y,  imax_x, imax_y, **select_opts2)
        )

    def update(self, start, end):
        # Current extrema of inner and outer rectangles.
        imin_x, imin_y,  imax_x, imax_y = self._get_coords(start, end)
        omin_x, omin_y,  omax_x, omax_y = 0, 0,  self.width, self.height

        # Update coords of all rectangles based on these extrema.
        self.canvas.coords(self.rects[0], omin_x, omin_y,  omax_x, imin_y),
        self.canvas.coords(self.rects[1], omin_x, imin_y,  imin_x, imax_y),
        self.canvas.coords(self.rects[2], imax_x, imin_y,  omax_x, imax_y),
        self.canvas.coords(self.rects[3], omin_x, imax_y,  omax_x, omax_y),
        self.canvas.coords(self.rects[4], imin_x, imin_y,  imax_x, imax_y),

        for rect in self.rects:  # Make sure all are now visible.
            self.canvas.itemconfigure(rect, state=tk.NORMAL)

    def _get_coords(self, start, end):
        """ Determine coords of a polygon defined by the start and
            end points one of the diagonals of a rectangular area.
        """
        return (min((start[0], end[0])), min((start[1], end[1])),
                max((start[0], end[0])), max((start[1], end[1])))

    def hide(self):
        for rect in self.rects:
            self.canvas.itemconfigure(rect, state=tk.HIDDEN)


class Application(tk.Frame):

    # Default selection object options.
    SELECT_OPTS = dict(dash=(2, 2), stipple='gray25', fill='red',
                          outline='')

    def __init__(self, parent, path, *args, **kwargs):
        super(parent, self).__init__(parent, *args, **kwargs)

        self.path = path
        img = ImageTk.PhotoImage(Image.fromarray(self.path))
        self.canvas = tk.Canvas(root, width=img.width(), height=img.height(), borderwidth=0, highlightthickness=0)
        self.canvas.pack(expand=True)

        self.canvas.create_image(0, 0, image=img, anchor=tk.NW)
        self.canvas.img = img  # Keep reference.

        # Create selection object to show current selection boundaries.
        self.selection_obj = SelectionObject(self.canvas, self.SELECT_OPTS)
        self.coords = []
        # Callback function to update it given two points of its diagonal.
        def on_drag(start, end, **kwarg):  # Must accept these arguments.
            self.selection_obj.update(start, end)
            self.coords.append(start)
            self.coords.append(end)

        # Create mouse position tracker that uses the function.
        self.posn_tracker = MousePositionTracker(self.canvas)
        self.posn_tracker.autodraw(command=on_drag)  # Enable callbacks.


# if __name__ == '__main__':
#
#     WIDTH, HEIGHT = 900, 900
#     BACKGROUND = 'grey'
#     TITLE = 'Image Cropper'
#
#     root = tk.Tk()
#     root.title(TITLE)
#     root.geometry('%sx%s' % (WIDTH, HEIGHT))
#     root.configure(background=BACKGROUND)
#     path = '/hpc/bsha219/campus-map.jpg'
#
#     app = Application(root, path, background=BACKGROUND)
#     app.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.TRUE)
#     app.mainloop()
#     print(app.coords[0][0])
#     print(app.coords[0][1])
#     print(app.coords[-1][0])
#     print(app.coords[-1][1])
#     print(app.posn_tracker.canv_width)
#     print(type(app.posn_tracker.canv_width))
#     # m = np.zeros((678, 577))
#     m = np.zeros((int(app.posn_tracker.canv_height), int(app.posn_tracker.canv_width)))
#     print(m.shape)
#     m[int(app.coords[0][1]):int(app.coords[-1][1]), int(app.coords[0][0]):int(app.coords[-1][0])] = 1
#     plt.imshow(m)
#     plt.show()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    # data_dir = '/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Raw/DICOMS'
    patient_scans = load_scan('/hpc/bsha219/_12Labour/Alfred/Imaging data/002/Post/scans/5-DE_CTPA_Fines__F_0.7/resources/DICOM/')
    patient_images = get_pixels_hu(patient_scans)
    filelist_mpa = glob.glob('/hpc/bsha219/_12Labour/Alfred/Imaging data/002/Post/scans/5-DE_CTPA_Fines__F_0.7/resources/DICOM/*')
    filelist_mpa.sort()
    global artery_mask
    artery_mask = []
    seed_found = False

    global coords
    coords = []
    WIDTH, HEIGHT = 900, 900
    BACKGROUND = 'grey'
    TITLE = 'Select ROI'

    # for k in range(len(filelist_mpa)):
    k = 0
    while seed_found == False:
        # for i in range(patient_images[k].shape[0]):
        #     for j in range(patient_images[k].shape[0]):
        #         if patient_images[k][i][j] > 0:
        root = tk.Tk()
        root.title(TITLE)
        root.geometry('%sx%s' % (WIDTH, HEIGHT))
        root.configure(background=BACKGROUND)

        app = Application(root, patient_images[k], background=BACKGROUND)
        app.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.TRUE)
        app.mainloop()
                    # patient_images[k][i][j] *= 12
    # blurred_f = ndimage.gaussian_filter(img, sigma=3)
        seed = np.zeros((len(patient_images), int(app.posn_tracker.canv_height), int(app.posn_tracker.canv_width)))
        if len(app.coords) > 0:
            seed[int(app.coords[0][1]):int(app.coords[-1][1]), int(app.coords[0][0]):int(app.coords[-1][0])] = 1
        filter_blurred_f = filters.median_filter(patient_images[k], 3)
        filter_blurred_f2 = ndimage.gaussian_filter(filter_blurred_f, sigma=1)
        alpha = 12
        patient_images[k] = filter_blurred_f + alpha * (filter_blurred_f - filter_blurred_f2)
        if seed.any() == 1:
            seed_found = True
        k = k + 1
    seg, phi, its = chanvese3d(patient_images, seed, max_its=1500, display=False, alpha=0.2)
    segment = seg * 5.0
    segment = np.uint8(segment)
    # to save them as a stack of masks in a directory
    binary_im = Image.fromarray(segment)
    binary_im.save('/hpc/bsha219/_12Labour/Alfred/Imaging data/002/Post/Vessel/TEST_MPA/TEST%.4d.jpg' % k,
                   quality=100)


