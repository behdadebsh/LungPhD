wss_file = '/hpc/bsha219/lung/Data/CTEPH/NewCase/FRC/CFD/3D/3/wallShearStress'
with open(wss_file) as f_in:
    lines = (line.rstrip() for line in f_in)
    wss_contents = list(line for line in lines if line)

num_points = int(wss_contents[24].split()[0])
# print(float(wss_contents[26].split()[0][1:-1]))
# x = float(wss_contents[28].split()[0][1:])
# y = float(wss_contents[28].split()[1])
# z = float(wss_contents[28].split()[2][0:-1])
# magnitude = (x ** 2 + y ** 2 + z ** 2) ** (0.5)
# print(wss_contents[27].split()[0][1:-1])
# print(x)
# print(z)
# print(magnitude)
magnitude = 0
for wss_line in range(num_points):
    x = float(wss_contents[26+wss_line].split()[0][1:])
    y = float(wss_contents[26+wss_line].split()[1])
    z = float(wss_contents[26+wss_line].split()[2][:-1])
    magnitude = magnitude + (x**2 + y**2 + z**2)**(0.5)
    # print(magnitude)

avg = magnitude/num_points
print(avg)

