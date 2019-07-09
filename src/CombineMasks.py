import matplotlib.image as mpimg
import numpy as np
from PIL import Image
import glob

filelist_mpa = glob.glob('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/MPA_Mask/*.png')
filelist_mpa.sort()
MPA_imgs = np.array([np.array(Image.open(fname)) for fname in filelist_mpa])
filelist_vessel = glob.glob('/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/InsideVessels/*.png')
filelist_vessel.sort()
vessel_imgs = np.array([np.array(Image.open(fname)) for fname in filelist_vessel])
GeneralMask = vessel_imgs + MPA_imgs
for slice in range(len(GeneralMask)):
    mpimg.imsave(
        '/hpc/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/whole_mask/FullMask%.4d.jpg' % slice,
        GeneralMask[slice])
