from load_scan import load_scan, get_pixels_hu
from Lung_Segmentation import lung_segment, generate_markers
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import matplotlib.image as im


ct_images = load_scan('../examples/sample_inputs')
images_hu = get_pixels_hu(ct_images)
for i in range(len(images_hu)):
    marker_internal, marker_watershed = generate_markers(images_hu[i])
    marker_internal[marker_internal == 2] = 1
    lungs, lung_left, lung_right = lung_segment(images_hu[i])
    lungs[lungs == 17] = 1
    lungs[lungs == 5] = 1
    ex = lungs - marker_internal*1.0
    vessel = ndimage.binary_opening(ex)
    plt.imshow(vessel)
    plt.show()
    im.imsave('./results/Vessel' + str(i+1) + '.jpg', vessel)
