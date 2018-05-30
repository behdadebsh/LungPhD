
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


