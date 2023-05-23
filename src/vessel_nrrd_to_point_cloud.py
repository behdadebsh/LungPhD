import numpy as np
import skimage.measure as measure
import nrrd
import os

def writeExDataFile(filename, coords, mean_radius_field=None):
    """
    Write out an ex data file with the given coords.
    :param filename: Filename to write to.
    :param coords: List of coordinate lists.
    :param mean_radius_field: Additonal mean radius field to write out (optional).
    :return: None
    """
    with open(filename, 'w') as f:
        f.write(' Group name : MPA_datapoints\n')
        f.write(' #Fields={0}\n'.format(1 if mean_radius_field is None else 2))
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('\tx.  Value index= 1, #Derivatives= 0\n')
        f.write('\ty.  Value index= 2, #Derivatives= 0\n')
        f.write('\tz.  Value index= 3, #Derivatives= 0\n')
        if mean_radius_field is not None:
            f.write(' 2) radius, field, rectangular cartesian, #Components=1\n')
            f.write('\tr.  Value index= 4, #Derivatives= 0\n')
        for i in range(len(coords)):
            f.write(' Node:     %.4d\n' % (i + 1))
            f.write('   %s  %s  %s\n' % (coords[i][0], coords[i][1], coords[i][2]))
            if mean_radius_field is not None:
                f.write('   %s\n' % mean_radius_field[i])

# Read the NRRD segmentation file
root_folder = '/hpc/bsha219/lung/Data'
study = 'CTEPH'
protocol = 'Post'
subject = 'Alfred1'
data, header = nrrd.read(os.path.join(root_folder, study, subject, protocol, 'Vessel/MPA_Mask', subject + '_' + protocol + '.nrrd'))
# data, header = nrrd.read('/hpc/bsha219/lung/Data/CTEPH/Alfred1/Pre/Alfred1_PRE_PRASH_FIJI_CENTERLINE.nrrd')
transposed = np.transpose(data, (2, 0, 1))
transposed = np.flip(transposed, axis=0)
transposed = np.transpose(transposed, (1, 2, 0))
xy_res = header['space directions'][0][0]
z_res = header['space directions'][2][2]

# Extract the external surface coordinates
verts, faces, _, _ = measure.marching_cubes(transposed)

# Down-sample the surface coordinates
num_samples = 3000
if len(verts) > num_samples:
    verts = verts[np.random.choice(len(verts), num_samples, replace=False)]

if z_res < 0:    # to make sure the orientation is consistent within all data
    verts = verts * [xy_res, xy_res, -z_res]
    verts = verts + [0, 0, data.shape[2]*z_res]
else:
    verts = verts * [xy_res, xy_res, -z_res]
# Write out the down-sampled coordinates to a file
writeExDataFile(os.path.join(root_folder, study, subject, protocol, 'Vessel/MPA_Mask', subject + '_' + protocol + '_points.exdata'), verts)

# set dir /hpc/bsha219/lung/Data/CTEPH/Alfred3/Post/Vessel/MPA_Mask
# gfx read data Alfred3_Post_points
