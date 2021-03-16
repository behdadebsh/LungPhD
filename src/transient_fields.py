import numpy as np


def readfieldfile(path, field):
    if field == 'flow' or 'pressure':
        fullpath = path + 'total_' + field + '.txt'
    if field == 'radius':
        fullpath = path + 'terminal_radii.txt'
    if field == 'WSS':
        fullpath = path + 'WSS.txt'
    # else:
    #     print("Field is not valid - you should choose between flow, pressure, radius and WSS")
    #     exit(0)
    f = open(fullpath, 'r')
    contents = f.readlines()
    l = []
    fields = np.zeros((len(contents),50),dtype=float)
    for line in range(len(contents)):
        l.append(contents[line].split())
    for line in range(len(l)):
        for timestep in range(50):
            fields[line][timestep] = l[line][timestep + 1]
    return fields


def splitlines(textlines):
    l = []
    for line in range(len(textlines)):
        l.append(textlines[line].split())
    return l


def readcoordinates(path, filename):
    """
    Reading the terminal nodes coordinates from file
    :param path: Path to directory with your terminal values exnode file
    :param filename: The filename you have saved your results as (usually the subject name)
    :return: coords: coordinates of terminal nodes which is a list of dictionaries (node_num,x,y,z)
    """
    f = open(path + filename + '.exnode', 'r')
    contents = f.readlines()
    coords = []
    for lines in range(len(contents)):
        line = contents[lines].split()
        if line[0] == 'Node:':
            coords.append(dict(node_num=line[1], x=contents[lines+1], y=contents[lines+2], z=contents[lines+3]))
    return coords


def writeExNodeFile(filename, coords, field, timestep):
    """
    Write out an ex Node file with the given coords.
    :param filename: Filename to write to.
    :param coords: List of coordinate lists.
    :return: an Exnode file
    """
    with open(filename, 'w') as f:
        f.write(' Group name : terminal radii\n')
        f.write(' #Fields=2\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('  x.  Value index=1, #Derivatives=0\n')
        f.write('  y.  Value index=2, #Derivatives=0\n')
        f.write('  z.  Value index=3, #Derivatives=0\n')
        f.write(' 2) flow, field, rectangular cartesian, #Components=1\n')
        f.write('  1.  Value index=4, #Derivatives=0\n')

        for i in range(len(coords)):
            f.write(' Node:           %d\n' % int(coords[i]['node_num']))
            f.write('     %s\n' % float(coords[i]['x']))
            f.write('     %s\n' % float(coords[i]['y']))
            f.write('     %s\n' % float(coords[i]['z']))
            f.write('     %s\n' % (field[i][timestep]))


def main():
    timestep = [0, 9, 19, 29, 39, 49]
    cycle_time = ['00', '16', '32', '48', '64', '80']
    flow_field = readfieldfile('/hpc/bsha219/hpcbuilds/CTEPH_terminal_comparison/CTEPH_20percent1_FlowBC/', 'flow')
    coords = readcoordinates('/hpc/bsha219/hpcbuilds/CTEPH_terminal_comparison/CTEPH_20percent1_FlowBC/','P2BRP268-H12816_terminal')
    i = 0
    for step in timestep:
        writeExNodeFile('/hpc/bsha219/hpcbuilds/CTEPH_terminal_comparison/CTEPH_20percent1_FlowBC/terminal_flow_' + cycle_time[i] + '.exnode', coords, flow_field, step)
        i = i + 1


if __name__ == '__main__':
    main()