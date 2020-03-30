# import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import copy
import math
# import scipy.integrate as integrate

# def U_prime(N):  # defining the gradient of potential(dimentionless for L=1. ):
#
#     U0 = 3
#     x = []
#     y = []
#     for i in range(N):
#         [a, b] = np.random.normal(0, 1, 2)
#         a = 6 * U0 * a ** 3 - a
#         b = 6 * U0 * b ** 3 - b
#         x.append(a)
#         y.append(b)
#     x = np.asarray(x)
#     y = np.asarray(y)
#
#     return x, y
#
#
# def main():
#     N = int(input("Please provide N:"))  # number of points for potential
#     NP = int(input("Please provide Number of particles:"))  # number of Particles involved
#     for j in range(NP):
#         x, y = U_prime(N)
#     print("x array is: ", x)
#     print("y array is: ", y)
#
#
# # if __name__ == '__main__':
# #     main()
# [a, b] = np.random.normal(0, 1, 2)
# print(a)
# print(b)
# dims = 2
# step_n = 10
# step_set = [-1, 0, 1]
# origin = np.zeros((1,dims))
# # Simulate steps in 2D
# step_shape = (step_n,dims)
# steps = np.random.choice(a=step_set, size=step_shape)
# path = np.concatenate([origin, steps]).cumsum(0)
# start = path[:1]
# stop = path[-1:]
# # Plot the path
# fig = plt.figure(figsize=(8,8),dpi=200)
# ax = fig.add_subplot(111)
# # ax.scatter(path[:,0], path[:,1],c='blue')
# ax.scatter(path[:,0], path[:,1],c='blue',alpha=0.25,s=0.05)
# ax.plot(path[:,0], path[:,1],c='blue',alpha=0.5,lw=0.25,ls='-')
# ax.plot(start[:,0], start[:,1],c='red', marker='+')
# ax.plot(stop[:,0], stop[:,1],c='black', marker='o')
# plt.title('2D Random Walk')
# plt.tight_layout(pad=0)
# plt.savefig('random_walk_2d.png',dpi=250)
