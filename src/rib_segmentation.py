import numpy as np
import pydicom as dicom
import os
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filters
import scipy.ndimage as ndimage
from skimage import segmentation, morphology
from mayavi import mlab
from skimage.morphology import cube, square
from Lung_Segmentation import lung_segment
from PIL import Image
from skimage.filters import frangi,hessian


def load_scan(path):
    slices = [dicom.read_file(path + '/' + s) for s in os.listdir(path)]
    slices.sort(key=lambda x: int(x.InstanceNumber))
    try:
        slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
    except:
        slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)
    for s in slices:
        s.SliceThickness = slice_thickness
    return slices


def get_pixels_hu(scans):
    """
    method to convert to hounsfield units.
    :param scans: slices from load_scan method.
    :return: a numpy array with hounsfield values
    """
    image = np.stack([s.pixel_array for s in scans])
    # Convert to int16 (from sometimes int16),
    # should be possible as values should always be low enough (<32k)
    image = image.astype(np.int16)
    # Set outside-of-scan pixels to 0
    # The intercept is usually -1024, so air is approximately 0
    image[image == -2000] = 0
    # Convert to Hounsfield units (HU)
    intercept = scans[0].RescaleIntercept
    slope = scans[0].RescaleSlope
    if slope != 1:
        image = slope * image.astype(np.float64)
        image = image.astype(np.int16)
    image += np.int16(intercept)
    return np.array(image, dtype=np.int16)


def generate_markers(image, threshold):
    # Creation of the internal Marker
    marker_internal = image < threshold  # Bone Tissue threshold
    marker_internal = segmentation.clear_border(marker_internal)
    return marker_internal


def cc_thresholding (mask):
    label_im, nb_labels = ndimage.label(mask)
    unique, counts, indices = np.unique(label_im, return_counts=True, return_index=True)
    # print unique, counts, indices
    threshold_labels = []
    for k in range(0, len(unique)):
        if indices[k] > 5:
            threshold_labels.append(unique[k])
    # print threshold_labels
    new_mask = np.where(np.isin(label_im, threshold_labels), label_im, 0)
    return new_mask


def main():
    path = '/hpc/bsha219/lung/Data/Human_Lung_Atlas/BRPILD001-H687/FRC/Raw'
    patient_scans = load_scan(path)
    patient_images = get_pixels_hu(patient_scans)
    # marker_internal = np.zeros([len(patient_scans),len(patient_images[0]), len(patient_images[0][0])])
    # bones = np.array([len(patient_images), 512, 512])
    bones = []
    segmented = []
    # for i in range(len(patient_images)): #this
    #     lungs, left_lung, right_lung = lung_segment(patient_images[i]) #this
    #     segmented.append(lungs) #this
    #     lungs = np.uint8(lungs) #this
        # to save LungMask as a stack of masks in a directory
        # binary_im = Image.fromarray(lungs) #this
        # binary_im.save('/hpc/bsha219/Hari/RA_work/BRPILD001-H687/FRC/Lung/LungMask%.4d.jpg' % i, quality=100)
        # filter_blurred = ndimage.gaussian_filter(patient_images[i], sigma=3)
        # marker_internal = generate_markers(patient_images[i], 100)  # you can play with threshold value to maybe get better results #this
        # marker_internal = morphology.binary_closing(marker_internal)#, selem=morphology.disk(1))
        # marker_internal = morphology.dilation(marker_internal)
        # marker_internal = marker_internal * 1
        # marker_internal = ndimage.binary_opening(marker_internal)
        # marker_internal = ndimage.binary_closing(marker_internal, structure=square(10)) #this
        # marker_internal = ndimage.binary_opening(marker_internal, structure=square(3))
        # marker_internal = ndimage.binary_opening(marker_internal)
        # marker_internal = ndimage.binary_closing(marker_internal)
        # marker_internal = marker_internal * 1 #this
        # print(marker_internal[30])
        # marker_internal = cc_thresholding(marker_internal)
        # marker_internal = morphology.opening(marker_internal)
        # bones.append(marker_internal) #this
    # to save them as a stack of masks in a directory
    #     plt.imshow(marker_internal)  # in case you want to see a slice
    #     plt.show()
    #     plt.imsave('/hpc/bsha219/Hari/RA_work/BRPILD001-H687/FRC/Ribcage_mask/BoneMask%.4d.png' % i, marker_internal)
    # bones = np.asarray(bones)
    # bones = np.flip(bones, axis=1)
    marker_internal = generate_markers(patient_images[100], 100)
    marker_internal = marker_internal * 1
    for i in range(len(patient_images)):
        ehem = frangi(patient_images[i], beta1=2, beta2=25, black_ridges=True)
    # ehem = frangi(patient_images[100])
    # ehem = hessian(patient_images[100])
    #     plt.imshow(ehem)
        plt.imsave('/hpc/bsha219/Hari/RA_work/BRPILD001-H687/FRC/Art_work/Frangi%.4d.png' % i, ehem)
    # plt.show()
    segmented = np.asarray(segmented)
    segmented = segmented * 1
    segmented = np.flip(segmented, axis=1)
    # print(type(bones))
    # print(bones.shape)
    # bones = ndimage.binary_opening(bones, structure=cube(3))
    # bones = bones * 1
    bones = ndimage.binary_closing(bones, structure=cube(10, dtype=np.uint8))
    bones = bones * 1
    # bones = ndimage.binary_opening(bones, structure=cube(3))
    # bones = bones * 1
    # bones = cc_thresholding(bones)
    bones = bones.T
    segmented = segmented.T
    # bones = morphology.binary_closing(bones, selem=cube(3, 3))
    # print(bones.shape)
    # print(type(bones))
    # mlab.mesh()
    # mlab.show()
    # src1 = mlab.pipeline.scalar_field(bones)
    # src2 = mlab.pipeline.scalar_field(segmented)
    # blur = mlab.pipeline.user_defined(src, filter='ImageGaussianSmooth')
    # voi = mlab.pipeline.extract_grid(blur)
    # mlab.pipeline.iso_surface(src1, colormap='Spectral')
    # mlab.pipeline.iso_surface(src2, colormap='Blues')
    # mlab.pipeline.volume(src1, vmin=0, vmax=0.8)
    # mlab.pipeline.image_plane_widget(src1, plane_orientation='x_axes',slice_index=10)
    # mlab.show()


if __name__ == '__main__':
    main()
