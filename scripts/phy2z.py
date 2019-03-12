import os
import sys
import numpy as np
import cantera as ct
import re

if len(sys.argv) != 3:
    print("Usage: python3 phy2z.py RawData.txt DstDir")
    exit(-1)
else:
    raw_data_path = sys.argv[1]
    output_dir = sys.argv[2]
    if not os.path.exists(raw_data_path):
        print("Input data file not found!")
        exit(-1)
    if not os.path.exists(output_dir):
        print("Output directory not found!")
        exit(-2)

gas = ct.Solution('gri30.cti')
MW = gas.molecular_weights
K = gas.n_species
P = ct.one_atm

CH4_IDX = gas.species_index('CH4')
H2_IDX = gas.species_index('H2')
O2_IDX = gas.species_index('O2')

def calc_rho(P, T, Y):
    ret = 0.0
    for k in range(K):
        ret += Y[k] / MW[k]
    ret = P / (ct.gas_constant * T * ret)
    return ret

def calc_mix_frac(Y):
    Yf = Y[CH4_IDX] + Y[H2_IDX]
    Yo = Y[O2_IDX]
    s = 40/9 # Mass ratio of oxidizer and fuel at stoichiometry
    Yo_0 = 0.232 # Mixture fraction of oxidizer in air stream
    Yf_0 = 1.0 # Mixture fraction of fuel in gas stream
    ret = (s*Yf-Yo+Yo_0)/(s*Yf_0+Yo_0) # The mixture fraction
    if ret < 0:
        ret = 0.0
    return ret

data = np.loadtxt(raw_data_path, skiprows=1)
N = len(data)
mdot_str = re.findall(r'\d+\.?\d*e?[-+]?\d+', raw_data_path)
mf = float(mdot_str[0])
mo = float(mdot_str[1])

x = data[:, 0]
u = data[:, 1]
V = data[:, 2]
T = data[:, 3]
A = data[0, 4]
Y = np.zeros([gas.n_species, N])
for k in range(gas.n_species):
    Y[k] = data[:, 5 + k]

Y_C = np.zeros(N)
Y_H = np.zeros(N)
Y_O = np.zeros(N)
Y_N = np.zeros(N)
MixFrac = np.zeros(N)
dZdn = np.zeros(N)

for n in range(N):
    loc_y = np.array([Y[k, n] for k in range(K)])
    gas.TPY = P, T[n], loc_y
    Y_C[n] = gas.elemental_mass_fraction('C')
    Y_H[n] = gas.elemental_mass_fraction('H')
    Y_O[n] = gas.elemental_mass_fraction('O')
    Y_N[n] = gas.elemental_mass_fraction('N')
    MixFrac[n] = calc_mix_frac(loc_y)


fn = "mf={:f}_mo={:f}_transformed.txt".format(mf, mo)
fp = os.path.join(output_dir, fn)

fout = open(fp, 'w')

header = "{:>18s}".format("rho")
header += "{:>18s}{:>18s}{:>18s}{:>18s}".format("Y_C", "Y_H", "Y_O", "Y_N")
header += "{:>18s}{:>18s}{:>18s}\n".format('Z', 'dZdn', 'kai')
fout.write(header)

for n in range(N):
    loc_y = np.array([Y[k, n] for k in range(K)])
    loc_density = calc_rho(P, T[n], loc_y)

    content = '{:>18.8e}'.format(loc_density)
    content += '{:>18.8e}{:>18.8e}{:>18.8e}{:>18.8e}'.format(Y_C[n], Y_H[n], Y_O[n], Y_N[n])
    content += '{:>18.8e}\n'.format(MixFrac[n])
    fout.write(content)
