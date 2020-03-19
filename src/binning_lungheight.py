

def readnodefile(path):
    f = open(path + '.exnode', 'r')
    contents = f.readlines()
    list = []
    for line in range(len(contents)):
        lines = contents[line].split()
        if lines[0] == 'Node:':
            list.append(dict(node_num=int(lines[1]), x=float(contents[line+1]), y=float(contents[line+2]), z=float(contents[line+3]), field=float(contents[line+4])))
    return list


def main():
    a = readnodefile('/hpc/bsha219/hpcbuilds/functional-models/wave_transmission_Ebrahimi2019/Remod_grade10/Remod_grade10/terminal_flow_80')
    direction = 'y'  # you can only choose the axis
    maximum = max(a, key=lambda x:x[direction])
    minimum = min(a, key=lambda x:x[direction])
    interval = abs(maximum[direction] - minimum[direction])/10
    a1 = minimum[direction] + interval
    a2 = a1 + interval
    a3 = a2 + interval
    a4 = a3 + interval
    a5 = a4 + interval
    a6 = a5 + interval
    a7 = a6 + interval
    a8 = a7 + interval
    a9 = a8 + interval

    cluster1 = []
    cluster2 = []
    cluster3 = []
    cluster4 = []
    cluster5 = []
    cluster6 = []
    cluster7 = []
    cluster8 = []
    cluster9 = []
    cluster10 = []
    for i in range(len(a)):
        if a[i][direction] < a1:
            cluster1.append(a[i])
        elif a1 < a[i][direction] < a2:
            cluster2.append(a[i])
        elif a2 < a[i][direction] < a3:
            cluster3.append(a[i])
        elif a3 < a[i][direction] < a4:
            cluster4.append(a[i])
        elif a4 < a[i][direction] < a5:
            cluster5.append(a[i])
        elif a5 < a[i][direction] < a6:
            cluster6.append(a[i])
        elif a6 < a[i][direction] < a7:
            cluster7.append(a[i])
        elif a7 < a[i][direction] < a8:
            cluster8.append(a[i])
        elif a8 < a[i][direction] < a9:
            cluster9.append(a[i])
        elif a9 < a[i][direction]:
            cluster10.append(a[i])

    print(float(sum(d['field'] for d in cluster1)) / len(cluster1))
    print(float(sum(d['field'] for d in cluster2)) / len(cluster2))
    print(float(sum(d['field'] for d in cluster3)) / len(cluster3))
    print(float(sum(d['field'] for d in cluster4)) / len(cluster4))
    print(float(sum(d['field'] for d in cluster5)) / len(cluster5))
    print(float(sum(d['field'] for d in cluster6)) / len(cluster6))
    print(float(sum(d['field'] for d in cluster7)) / len(cluster7))
    print(float(sum(d['field'] for d in cluster8)) / len(cluster8))
    print(float(sum(d['field'] for d in cluster9)) / len(cluster9))
    print(float(sum(d['field'] for d in cluster10)) / len(cluster10))

if __name__ == '__main__':
    main()