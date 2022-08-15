#!/usr/bin/python
import numpy as np

class FlowIntensity(object):

    def __init__(self, dicom_dir, mask_dir):
        import glob
        from PIL import Image
        import pydicom as dicom
        import os

        self.dicoms = []
        self.masks = []
        self.filelist_dicoms = glob.glob(dicom_dir + '*')
        self.filelist_dicoms.sort()
        self.filelist_masks = glob.glob(mask_dir + '*.jpg')
        self.filelist_masks.sort()
        # load the DICOM files
        # print('glob: {}'.format(sys.argv[1]))
        self.slices = [dicom.read_file(dicom_dir + '/' + s) for s in os.listdir(dicom_dir)]
        self.slices.sort(key=lambda x: int(x.InstanceNumber))
        try:
            slice_thickness = np.abs(self.slices[0].ImagePositionPatient[2] - self.slices[1].ImagePositionPatient[2])
        except:
            slice_thickness = np.abs(self.slices[0].SliceLocation - self.slices[1].SliceLocation)

        for s in self.slices:
            s.SliceThickness = slice_thickness

        self.image = np.stack([s.pixel_array for s in self.slices])
        # Convert to int16 (from sometimes int16),
        # should be possible as values should always be low enough (<32k)
        self.image = self.image.astype(np.int16)

        # Set outside-of-scan pixels to 0
        # The intercept is usually -1024, so air is approximately 0
        self.image[self.image == -2000] = 0

        # Convert to Hounsfield units (HU)
        intercept = self.slices[0].RescaleIntercept
        slope = self.slices[0].RescaleSlope

        if slope != 1:
            image = slope * self.image.astype(np.float64)
            image = image.astype(np.int16)

        self.image += np.int16(intercept)
        # load masks
        for fname in self.filelist_masks:
            self.masks.append(np.array(Image.open(fname)))


mask = '/hpc/bsha219/lung/Data/Human_Lung_Atlas/P2BRP268-H12816/FRC/Vessel/'
dicoms = '/hpc/bsha219/lung/Data/Human_Lung_Atlas/P2BRP268-H12816/FRC/DICOMS/'
a = FlowIntensity(dicoms, mask)
