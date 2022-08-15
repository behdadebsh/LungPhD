import numpy as np
import matplotlib.pyplot as plt
from numpy import cos, sin, sqrt, pi

path = '/hpc/bsha219/Packages/functional-models/wave_transmission_Ebrahimi2019/Remod_grade1_FlowBC83338/'
# path1 = '/hpc/bsha219/Packages/functional-models/wave_transmission_Ebrahimi2019/Remod_grade1/'
flow_path = path + 'Inlet_flow.txt'
pressure_path = path + 'Inlet_pressure.txt'
forward_flow_path = path + 'Inlet_forward_flow.txt'
reflected_flow_path = path + 'Inlet_reflected_flow.txt'
forward_press_path = path + 'Inlet_forward_pressure.txt'
reflected_press_path = path + 'Inlet_reflected_pressure.txt'
flow_file = open(flow_path, 'r')
pressure_file = open(pressure_path, 'r')
forward_flow_file = open(forward_flow_path, 'r')
reflected_flow_file = open(reflected_flow_path, 'r')
forward_press_file = open(forward_press_path, 'r')
reflected_press_file = open(reflected_press_path, 'r')
flows = flow_file.readlines()
pressures = pressure_file.readlines()
forward_flows = forward_flow_file.readlines()
reflected_flows = reflected_flow_file.readlines()
forward_press = forward_press_file.readlines()
reflected_press = reflected_press_file.readlines()
flow_timesteps = flows[0].split()
pressure_timesteps = pressures[0].split()
forward_flow_timesteps = forward_flows[0].split()
reflected_flow_timesteps = reflected_flows[0].split()
forward_press_timesteps = forward_press[0].split()
reflected_press_timesteps = reflected_press[0].split()
flow = np.zeros(50, dtype=float)
pressure = np.zeros(50, dtype=float)
f_flow = np.zeros(50, dtype=float)
r_flow = np.zeros(50, dtype=float)
f_press = np.zeros(50, dtype=float)
r_press = np.zeros(50, dtype=float)
for i in range(50):
    flow[i] = float(flow_timesteps[i])
    pressure[i] = float(pressure_timesteps[i])
    f_flow[i] = float(forward_flow_timesteps[i])
    r_flow[i] = float(reflected_flow_timesteps[i])
    f_press[i] = float(forward_press_timesteps[i])
    r_press[i] = float(reflected_press_timesteps[i])


A0 = 2.8795/2
t = np.linspace(0.0, 0.8, 50)
w = 2 * pi * 72/60
# A = np.array([A0, 1.8971, 1.3188, 0.3636, 0.3037, 0.0541, 0.1061, 0.1004, 0.0121, 0.0314, 0.0197])
A = np.array([10.4244, 5.4375, 3.0769, 0.5502, 0.5689, 0.8515, 0.4410, 0.2950])
# A = A*100000
# phi = np.array([-1.8550, 2.5083, 0.204, 3.1331, 0.9108, 2.7929, 1.2929, 1.9993, 2.4006, -2.1450])
phi = np.array([-1.9528, -2.5565, 2.9375, 1.8215, -1.3285, -2.4211, 2.9831, -1.1246])
q = A[0]*cos(w*t + phi[0]) + A[1]*cos(2*w*t + phi[1]) + A[2]*cos(3*w*t + phi[2]) + A[3]*cos(4*w*t + phi[3]) + A[4]*cos(5*w*t + phi[4]) + A[5]*cos(6*w*t + phi[5]) + A[6]*cos(7*w*t + phi[6]) + A[7]*cos(8*w*t + phi[7])# + A[9]*cos(9*w*t + phi[8]) + A[10]*cos(20*w*t + phi[9])
# q = q * 6 / 100000
time = np.arange(0.0, 0.8, 0.016)
# fig, ax = plt.subplots(2, 1)
plt.figure()
plt.plot(t, q)
plt.show()
pressure = pressure * 0.00750062  # Pa to mmHg
f_press = f_press * 0.00750062  # Pa to mmHg
r_press = r_press * 0.00750062  # Pa to mmHg
flow = flow * 6 * 10**(-5)  # mm3/s to l/min
f_flow = f_flow * 6 * 10**(-5)  # mm3/s to l/min
r_flow = r_flow * 6 * 10**(-5)  # mm3/s to l/min
plt.figure()
plt.plot(time, flow, color='r')
plt.xlabel('Time [sec]')
plt.ylabel('Inlet Total flow '+r'$ ( \frac{l}{min}$)')
plt.show()
ax1 = plt.subplot(321)
ax1.set_ylabel('Inlet Total Pressure [mmHg]')
# ax1.set_ylim(-5, 30)
plt.plot(time, pressure, color='tab:blue', marker='o')
ax2 = plt.subplot(322)
ax2.set_ylabel('Inlet Total flow '+r'$ ( \frac{l}{min}$)')
# ax2.set_xlabel('Time [sec]')
plt.plot(time, flow, color='r')
ax5 = plt.subplot(323)
ax5.set_ylabel('Forward Pressure [mmHg]')
# ax5.set_ylim(-5, 30)
# ax5.set_xlabel('Time [sec]')
plt.plot(time, f_press, color='tab:blue', marker='o')
ax3 = plt.subplot(324)
ax3.set_ylabel('Forward Flow '+r'$ ( \frac{l}{min}$)')
# ax3.set_xlabel('Time [sec]')
plt.plot(time, f_flow, color='r')
ax6 = plt.subplot(325)
ax6.set_ylabel('Reflected Pressure [mmHg]')
# ax6.set_ylim(-5, 30)
ax6.set_xlabel('Time [sec]')
plt.plot(time, r_press, color='tab:blue', marker='o')
ax4 = plt.subplot(326)
ax4.set_ylabel('Reflected Flow '+r'$ ( \frac{l}{min}$)')
ax4.set_xlabel('Time [sec]')
plt.plot(time, r_flow, color='r')
plt.show()
