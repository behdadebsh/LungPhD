import numpy as np

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
    # a = readnodefile('/hpc/bsha219/hpcbuilds/functional-models/wave_transmission_Ebrahimi2019/Remod_grade10/Remod_grade10/terminal_flow_80')
    a = readnodefile('E:\\wave_transmission\\Remod_grade10\\Remod_grade10\\terminal_flow_64')
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

    val1 = []
    val2 = []
    val3 = []
    val4 = []
    val5 = []
    val6 = []
    val7 = []
    val8 = []
    val9 = []
    val10 = []
    for in_cluster in range(len(cluster1)):
        val1.append(cluster1[in_cluster]['field'])
    for in_cluster in range(len(cluster2)):
        val2.append(cluster2[in_cluster]['field'])
    for in_cluster in range(len(cluster3)):
        val3.append(cluster3[in_cluster]['field'])
    for in_cluster in range(len(cluster4)):
        val4.append(cluster4[in_cluster]['field'])
    for in_cluster in range(len(cluster5)):
        val5.append(cluster5[in_cluster]['field'])
    for in_cluster in range(len(cluster6)):
        val6.append(cluster6[in_cluster]['field'])
    for in_cluster in range(len(cluster7)):
        val7.append(cluster7[in_cluster]['field'])
    for in_cluster in range(len(cluster8)):
        val8.append(cluster8[in_cluster]['field'])
    for in_cluster in range(len(cluster9)):
        val9.append(cluster9[in_cluster]['field'])
    for in_cluster in range(len(cluster10)):
        val10.append(cluster10[in_cluster]['field'])
    std1 = np.std(np.asarray(val1))
    std2 = np.std(np.asarray(val2))
    std3 = np.std(np.asarray(val3))
    std4 = np.std(np.asarray(val4))
    std5 = np.std(np.asarray(val5))
    std6 = np.std(np.asarray(val6))
    std7 = np.std(np.asarray(val7))
    std8 = np.std(np.asarray(val8))
    std9 = np.std(np.asarray(val9))
    std10 = np.std(np.asarray(val10))
    print('Standarad deviation for cluster 1:', std1)
    print('Average flow for cluster 1:', float(sum(d['field'] for d in cluster1)) / len(cluster1))
    print('Standarad deviation for cluster 2:', std2)
    print('Average flow for cluster 2:', float(sum(d['field'] for d in cluster2)) / len(cluster2))
    print('Standarad deviation for cluster 3:', std3)
    print('Average flow for cluster 3:', float(sum(d['field'] for d in cluster3)) / len(cluster3))
    print('Standarad deviation for cluster 4:', std4)
    print('Average flow for cluster 4:', float(sum(d['field'] for d in cluster4)) / len(cluster4))
    print('Standarad deviation for cluster 5:', std5)
    print('Average flow for cluster 5:', float(sum(d['field'] for d in cluster5)) / len(cluster5))
    print('Standarad deviation for cluster 6:', std6)
    print('Average flow for cluster 6:', float(sum(d['field'] for d in cluster6)) / len(cluster6))
    print('Standarad deviation for cluster 7:', std7)
    print('Average flow for cluster 7:', float(sum(d['field'] for d in cluster7)) / len(cluster7))
    print('Standarad deviation for cluster 8:', std8)
    print('Average flow for cluster 8:', float(sum(d['field'] for d in cluster8)) / len(cluster8))
    print('Standarad deviation for cluster 9:', std9)
    print('Average flow for cluster 9:', float(sum(d['field'] for d in cluster9)) / len(cluster9))
    print('Standarad deviation for cluster 10:', std10)
    print('Average flow for cluster 10:', float(sum(d['field'] for d in cluster10)) / len(cluster10))

if __name__ == '__main__':
    main()