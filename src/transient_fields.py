import numpy as np


def readfieldfile(path, field):
    if field != 'flow':
        field = 'pressure'
    fullpath = path + 'total_' + field + '.txt'
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
    :return: None
    """
    with open(filename, 'w') as f:
        f.write(' Group name : flow\n')
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
    flow_field = readfieldfile('/hpc/bsha219/hpcbuilds/functional-models/wave_transmission_Ebrahimi2019/Remod_grade10/Remod_grade10/', 'flow')
    coords = readcoordinates('/hpc/bsha219/hpcbuilds/functional-models/wave_transmission_Ebrahimi2019/Remod_grade10/Remod_grade10/','P2BRP268-H12816_terminal')
    i = 0
    for step in timestep:
        writeExNodeFile('/hpc/bsha219/hpcbuilds/functional-models/wave_transmission_Ebrahimi2019/Remod_grade10/Remod_grade10/terminal_flow_' + cycle_time[i] + '.exnode', coords, flow_field, step)
        i = i + 1



if __name__ == '__main__':
    main()