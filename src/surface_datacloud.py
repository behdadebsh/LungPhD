"""

Read in image masks and extracting surface data cloud

@author: Behdad Shaarbaf Ebrahimi (UoA - ABI) for academic purposes.
"""
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import glob
from utils import writeExDataFile

filelist_mpa = glob.glob('/hpc/bsha219/lung/Data/Human_Lung_Atlas/P2BRP268-H12816/FRC/Vessel/*.jpg')
filelist_mpa.sort()
MPA_imgs = np.array([np.array(Image.open(fname)) for fname in filelist_mpa])
coords = []
for slice_number in range(len(filelist_mpa)):
    image = plt.imread(filelist_mpa[slice_number])
    for y in range(image.shape[0]-1):
        FOUND = False
        for x in range(image.shape[1]-1):
            if not FOUND and (int(image[x+1][y])) == 1 and (int(image[x][y])) == 0:
                FOUND = True
                coords.append([(x+1)*0.65, y*0.65, slice_number*1.25])
            elif FOUND and (int(image[x+1][y])) == 0 and (int(image[x][y])) == 1:
                FOUND = False
                coords.append([x*0.65, y*0.65, slice_number*1.25])
writeExDataFile('/hpc/bsha219/lung/Data/Human_Lung_Atlas/P2BRP268-H12816/Vessel/surface_MPAtrimmed.exdata', coords, mean_radius_field=None)
