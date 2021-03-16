
def writeExDataFile(filename, coords, mean_radius_field=None):
    """
    Write out an ex data file with the given coords.
    :param filename: Filename to write to.
    :param coords: List of coordinate lists.
    :param mean_radius_field: Additonal mean radius field to write out (optional).
    :return: None
    """
    with open(filename, 'w') as f:
        f.write(' Group name : bif_centerline_points\n')
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


def writeExNodeFile(filename, coords, fields, list_of_fields):
    """
    Write out an ex Node file with the given coords.
    :param filename: Filename to write to.
    :param coords: List of coordinate lists.
    :return: None
    """
    with open(filename, 'w') as f:
        f.write(' Group name : MAC\n')
        f.write(' #Fields=%1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 0\n')
        f.write('   y.  Value index= 2, #Derivatives= 0\n')
        f.write('   z.  Value index= 3, #Derivatives= 0\n')
        for i in range(len(coords)):
            f.write(' Node:         %.4d\n' % (i + 1001))
            f.write('    %s\n' % (coords[i][0]))
            f.write('    %s\n' % (coords[i][1]))
            f.write('    %s\n' % (coords[i][2]))


def writeipNodeFile(filename, coords):
    """
    Write out an ipnode file with the given coords.
    :param filename: Filename to write to.
    :param coords: List of coordinate lists.
    :return: None
    """
    with open(filename, 'w') as f:
        f.write(' CMISS Version 2.1  ipnode File Version 2\n')
        f.write(' Heading: MAC\n\n')
        f.write(' The number of nodes is [    69]:     69\n')
        f.write(' Number of coordinates [3]: 3\n')
        f.write(' Do you want prompting for different versions of nj=1 [N]? N\n')
        f.write(' Do you want prompting for different versions of nj=2 [N]? N\n')
        f.write(' Do you want prompting for different versions of nj=3 [N]? N\n')
        f.write(' The number of derivatives for coordinate 1 is [0]: 0\n')
        f.write(' The number of derivatives for coordinate 2 is [0]: 0\n')
        f.write(' The number of derivatives for coordinate 3 is [0]: 0\n\n')
        for i in range(len(coords)):
            f.write(' Node number [  %d]:   %d\n' % ((i + 1001), (i + 1001)))
            f.write(' The Xj(1) coordinate is [ 0.00000E+00]:  %s\n' % (coords[i][0]))
            f.write(' The Xj(2) coordinate is [ 0.00000E+00]:  %s\n' % (coords[i][1]))
            f.write(' The Xj(3) coordinate is [ 0.00000E+00]:  %s\n\n' % (coords[i][2]))
