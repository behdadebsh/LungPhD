rm1_sol_path = '/hpc/bsha219/Research related/EMBC2023/Perfusion/RM1/block_mesh_flow_perf.exelem'
rm7_sol_path = '/hpc/bsha219/Research related/EMBC2023/Perfusion/RM7/block_mesh_flow_perf.exelem'

rm7_flows = []
rm1_flows = []
flow_fraction = []

with open(rm1_sol_path) as f_in:
    lines = (line.rstrip() for line in f_in)
    rm1_contents = list(line for line in lines if line)
with open(rm7_sol_path) as f_in:
    lines = (line.rstrip() for line in f_in)
    rm7_contents = list(line for line in lines if line)
for line in range(len(rm1_contents)):
    if rm1_contents[line].split()[0] == 'Element:':
        rm1_flows.append(float(rm1_contents[line + 2].split()[0]))
        rm7_flows.append(float(rm7_contents[line + 2].split()[0]))
        flow_fraction.append((float(rm7_contents[line + 2].split()[0]) - float(rm1_contents[line + 2].split()[0]))/float(rm1_contents[line + 2].split()[0]))

with open("/hpc/bsha219/Research related/EMBC2023/Intensity_mapping/flow_diff_fraction.exelem", 'w') as f:   # write to file in order: #Data_point,  #Terminal_node, X, Y, Z, Cluster
    f.write(" Group name: flow_diff1\n Shape.  Dimension=1\n #Scale factor sets=0\n #Nodes=0\n #Fields=1\n")
    f.write(" 1) flow_diff_fraction_RM1, field, rectangular cartesian, #Components=1\n")
    f.write("   flow_diff_fraction_RM1.  l.Lagrange, no modify, grid based.\n   #xi1=1\n")
    for i in range(len(flow_fraction)):
        f.write(" Element:     %d 0 0\n" % (i+1))
        f.write("   Values:\n")
        f.write("       %s    %s\n" % (flow_fraction[i], flow_fraction[i]))
        # if (term_field[i][5] == 7):
        #     f.write(' %d   %s   %s   %s  1.0  1.0  1.0\n' % (int(term_field[i][0]), float(term_field[i][2]),
        #                                           float(term_field[i][3]), float(term_field[i][4])))
f.close()
