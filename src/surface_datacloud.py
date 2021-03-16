"""

Read in image masks and extracting surface data cloud

@author: Behdad Shaarbaf Ebrahimi (UoA - ABI) for academic purposes.
"""
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import glob
from utils import writeExDataFile

# filelist_mpa = glob.glob('E:\lung\Data\Human_PE_Study_HRC\ST12\TLC\Vessel\MPA_Mask\*.png')
filelist_mpa = glob.glob('E:\lung\Data\CTEPH\CTEPH6\FRC\Vessel\MPA2\*.jpg')
filelist_mpa.sort()
MPA_imgs = np.array([np.array(Image.open(fname)) for fname in filelist_mpa])
# for slice in range(len(MPA_imgs)):
#     mpimg.imsave(
#         '/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/whole_mask/FullMask%.4d.jpg' % slice,
#         MPA_imgs[slice])
# image = plt.imread(filelist_mpa[1])
# print(image.shape[0])
coords = []
for slice_number in range(len(filelist_mpa)):
    image = plt.imread(filelist_mpa[slice_number])
    for y in range(image.shape[0]-1):
        FOUND = False
        for x in range(image.shape[1]-1):
            if not FOUND and (int(image[x+1][y][0])) == 0 and (int(image[x][y][0])) == 1:
                FOUND = True
                coords.append([x+1, y, slice_number])
            elif FOUND and (int(image[x+1][y][0])) == 1 and (int(image[x][y][0])) == 0:
                FOUND = False
                coords.append([x, y, slice_number])
writeExDataFile('E:\lung\Data\CTEPH\CTEPH6\FRC\Vessel\surface_MPAtrimmed.exdata', coords, mean_radius_field=None)
