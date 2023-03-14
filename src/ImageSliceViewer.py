# import numpy as np
# import matplotlib.pyplot as plt


# Fixing random state for reproducibility
# np.random.seed(19680801)
import pydicom
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider

# Load the DICOM file
ds = pydicom.dcmread("CT.dcm")

# Get the pixel dataPredicting a patient's response to pulmonary thromboendarterectomy (PTE) surgery for CTEPH is essential in clinical decision-making, as it helps to determine the best course of treatment and ensures the best possible outcomes. This can include assessing factors such as the extent of clotting, lung function, and overall health, and can aid in determining the likelihood of successful surgical removal of the clots. Accurately predicting response to PTE surgery can lead to improved patient outcomes and better quality of life for individuals with CTEPH.
pixel_data = ds.pixel_array

# Create a figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set the axis labels
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# Plot the pixel data as a stack of 3 perpendicular planes
image = ax.imshow(pixel_data[:,:,0], cmap=plt.cm.gray, extent=[0, ds.Rows, 0, ds.Columns, 0, ds.NumberOfFrames])

# Enable interactive mouse rotation and zoom
plt.subplots_adjust(left=0.25, bottom=0.25)
for i in range(3):
    ax.mouse_init(rotate_btn=i+1, zoom_btn=4)

# Add a slider to scroll through the planes
axcolor = 'lightgoldenrodyellow'
axslice = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
slice = Slider(axslice, 'Slice', 0, pixel_data.shape[2]-1, valinit=0, valstep=1)

def update(val):
    z = int(slice.val)
    image.set_data(pixel_data[:,:,z])
    ax.set_zlim(z, z+1)
    fig.canvas.draw_idle()
slice.on_changed(update)

plt.show()

# class IndexTracker:
#     def __init__(self, ax, X):
#         self.ax = ax
#         # ax.set_title('use scroll wheel to navigate images')
#
#         self.X = X
#         self.slices, rows, cols = X.shape
#         self.ind = self.slices//2
#
#         self.im = ax.imshow(self.X[self.ind, :, :])
#         self.update()
#
#     def on_scroll(self, event):
#         print("%s %s" % (event.button, event.step))
#         if event.button == 'up':
#             self.ind = (self.ind + 1) % self.slices
#         else:
#             self.ind = (self.ind - 1) % self.slices
#         self.update()
#
#     def update(self):
#         self.im.set_data(self.X[self.ind, :, :])
#         self.ax.set_ylabel('slice %s' % self.ind)
#         self.im.axes.figure.canvas.draw()


# fig, ax = plt.subplots(1, 1)
#
# X = np.random.rand(20, 20, 40)
#
# tracker = IndexTracker(ax, X)
#
#
# fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
# plt.show()