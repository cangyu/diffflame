import os
import numpy as np
from scipy import integrate
import cantera as ct
from matplotlib import pyplot as plt
import re

data_file_name = 'mf=1.800000e-02_mo=5.400000e-02_stat.txt'
data_file_path = os.path.join('.', data_file_name)
data = np.loadtxt(data_file_path, skiprows=1)
N, _ = data.shape

mdot_str = re.findall(r'\d+\.?\d*e?[-+]?\d+', data_file_name)
mf = float(mdot_str[0])
mo = -float(mdot_str[1])

z = data[:, 0]
rho = data[:, 1]
u = data[:, 2]
V = data[:, 3]
T = data[:, 4]
Y_N2 = data[:, -10]
Y_CH4 = data[:, 18]
Y_CO = data[:, 19]
Y_CO2 = data[:, 20]

Y_C = data[:, -4]
Y_H = data[:, -3]
Y_O = data[:, -2]
Y_N = data[:, -1]

#plt.plot(z, Y_C, label='Y_C')
#plt.plot(z, Y_H, label='Y_H')
#plt.plot(z, Y_O, label='Y_O')
#plt.plot(z, Y_N, label='Y_N')
#plt.plot(z, Y_CO, label='Y_CO')
plt.plot(z, Y_CH4, label='Y_CH4')
#plt.plot(z, Y_CO2, label='Y_CO2')
plt.ylabel('Y')
plt.legend()
plt.twinx()
plt.plot(z, T, 'r')
plt.ylabel('T/K')
plt.show()

starting_idx = 500
ending_idx = 3251
print('Starting position: z = {}m'.format(z[starting_idx]))
print('Ending position: z = {}m'.format(z[ending_idx-1]))
mdot_left = rho[starting_idx] * u[starting_idx]
mdot_right = rho[ending_idx] * u[ending_idx]

print('Checking total mass...')
LHS = mdot_left - mdot_right
print(LHS)
integrand = np.array([2*rho[i]*V[i] for i in range(starting_idx, ending_idx)])
znew = np.array([z[i] for i in range(starting_idx, ending_idx)])
RHS = integrate.simps(integrand, znew)
print(RHS)

print('Checking N2...')
LHS = mdot_left * 0.0 - mdot_right * 0.7647
print(LHS)
integrand = np.array([2*rho[i]*Y_N2[i]*V[i] for i in range(starting_idx, ending_idx)])
znew = np.array([z[i] for i in range(starting_idx, ending_idx)])
RHS = integrate.simps(integrand, znew)
print(RHS)

print('Checking C...')
LHS = mdot_left * 0.6667 - mdot_right * 0.0
print(LHS)
integrand = np.array([2*rho[i]*Y_C[i]*V[i] for i in range(starting_idx, ending_idx)])
znew = np.array([z[i] for i in range(starting_idx, ending_idx)])
RHS = integrate.simps(integrand, znew)
print(RHS)

print('Checking H...')
LHS = mdot_left * 0.3333 - mdot_right * 0.0
print(LHS)
integrand = np.array([2*rho[i]*Y_H[i]*V[i] for i in range(starting_idx, ending_idx)])
znew = np.array([z[i] for i in range(starting_idx, ending_idx)])
RHS = integrate.simps(integrand, znew)
print(RHS)

print('Checking O...')
LHS = mdot_left * 0.0 - mdot_right * 0.2331
print(LHS)
integrand = np.array([2*rho[i]*Y_O[i]*V[i] for i in range(starting_idx, ending_idx)])
znew = np.array([z[i] for i in range(starting_idx, ending_idx)])
RHS = integrate.simps(integrand, znew)
print(RHS)
