import os
import sys
import numpy as np
import cantera as ct


def relaxation(a, b, x: float):
    return (1 - x) * a + x * b


def df_central(f, x):
    pnt_num = len(x)
    ret = np.zeros(pnt_num)
    ret[0] = ((x[1] - x[0]) / (x[2] - x[0]) * (f[2] - f[0]) - (x[2] - x[0]) / (x[1] - x[0]) * (f[1] - f[0])) / (
            x[1] - x[2])
    for i in range(1, pnt_num - 1):
        ret[i] = ((x[i] - x[i - 1]) / (x[i + 1] - x[i]) * (f[i + 1] - f[i]) - (x[i + 1] - x[i]) / (x[i] - x[i - 1]) * (
                f[i - 1] - f[i])) / (x[i + 1] - x[i - 1])
    ret[-1] = ((x[-1] - x[-3]) / (x[-1] - x[-2]) * (f[-2] - f[-1]) - (x[-1] - x[-2]) / (x[-1] - x[-3]) * (
            f[-3] - f[-1])) / (x[-3] - x[-2])
    return ret


def df_upwind(f, x, upwind_indicator):
    pnt_num = len(x)
    ret = np.zeros(pnt_num)
    ret[0] = (f[1] - f[0]) / (x[1] - x[0])
    for i in range(1, pnt_num - 1):
        s = -1 if upwind_indicator[i] > 0 else 1
        ret[i] = (f[i + s] - f[i]) / (x[i + s] - x[i])
    ret[-1] = (f[-1] - f[-2]) / (x[-1] - x[-2])
    return ret


def ddf(f, x):
    pnt_num = len(x)
    ret = np.zeros(pnt_num)
    ret[0] = 2.0 / (x[2] - x[1]) * ((f[2] - f[0]) / (x[2] - x[0]) - (f[1] - f[0]) / (x[1] - x[0]))
    for i in range(1, pnt_num - 1):
        ret[i] = 2.0 / (x[i + 1] - x[i - 1]) * (
                (f[i - 1] - f[i]) / (x[i] - x[i - 1]) + (f[i + 1] - f[i]) / (x[i + 1] - x[i]))
    ret[-1] = 2.0 / (x[-3] - x[-2]) * ((f[-2] - f[-1]) / (x[-1] - x[-2]) - (f[-3] - f[-1]) / (x[-1] - x[-3]))
    return ret


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
gas.transport_model = 'UnityLewis'
MW = gas.molecular_weights
K = gas.n_species
P = ct.one_atm

CH4_IDX = gas.species_index('CH4')
H2_IDX = gas.species_index('H2')
O2_IDX = gas.species_index('O2')
N2_IDX = gas.species_index('N2')
AR_IDX = gas.species_index('AR')

MW_C = gas.atomic_weight('C')
MW_H = gas.atomic_weight('H')
MW_O = gas.atomic_weight('O')


def calc_density(P, T, Y):
    ret = 0.0
    for k in range(K):
        ret += Y[k] / MW[k]
    ret = P / (ct.gas_constant * T * ret)
    return ret



'''
Gas: 1 CH4 + 1 H2
Air: 0.78 N2 + 0.21 O2 + 0.01 AR
'''
# Mass ratio of oxidizer and fuel at stoichiometry
s = 40/9 
# Mass fraction of fuel in gas stream
Yf_fu = 1.0 
# Mass fraction of fuel in air stream
Yf_ox = 0.0 
# Mass fraction of C in gas stream
Yc_fu = MW_C / (MW[CH4_IDX] + MW[H2_IDX]) 
# Mass fraction of C in air stream
Yc_ox = 0.0 
# Mass fraction of H in gas stream
Yh_fu = 6 * MW_H / (MW[CH4_IDX] + MW[H2_IDX]) 
# Mass fraction of H in air stream
Yh_ox = 0.0 
# Mass fraction of O in gas stream
Yo_fu = 0.0 
# Mass fraction of O in air stream
Yo_ox = 0.21 * MW[O2_IDX] / (0.21 * MW[O2_IDX] + 0.78 * MW[N2_IDX] + 0.01 * MW[AR_IDX]) 


def calc_mix_frac(Y):
    Yf = Y[CH4_IDX] + Y[H2_IDX]
    Yo = Y[O2_IDX]
    ret = (s*Yf-Yo+Yo_ox)/(s*Yf_fu+Yo_ox) # The mixture fraction
    return max(ret, 0.0)


def bilger(Yc, Yh, Yo):
    a = 2 * (Yc - Yc_ox) / MW_C + (Yh - Yh_ox) / (2*MW_H) - 2 * (Yo - Yo_ox) / MW_O
    b = 2 * (Yc_fu - Yc_ox) / MW_C + (Yh_fu - Yh_ox) / (2*MW_H) - 2 * (Yo_fu - Yo_ox) / MW_O
    ret = a / b
    return max(ret, 0.0)


data = np.loadtxt(raw_data_path, skiprows=1)
N = len(data)

x = data[:, 0]
u = data[:, 1]
V = data[:, 2]
T = data[:, 3]
A = data[0, 4]
Y = np.zeros([gas.n_species, N])
for k in range(gas.n_species):
    Y[k] = data[:, 5 + k]

rho = np.zeros(N)
Y_C = np.zeros(N)
Y_H = np.zeros(N)
Y_O = np.zeros(N)
Y_N = np.zeros(N)
MixFrac = np.zeros(N)
dZdn = np.zeros(N)
D = np.zeros(N)
for n in range(N):
    loc_y = np.array([Y[k, n] for k in range(K)])
    gas.TPY = T[n], P, loc_y
    rho[n] = calc_density(P, T[n], loc_y)
    Y_C[n] = gas.elemental_mass_fraction('C')
    Y_H[n] = gas.elemental_mass_fraction('H')
    Y_O[n] = gas.elemental_mass_fraction('O')
    Y_N[n] = gas.elemental_mass_fraction('N')
    # MixFrac[n] = calc_mix_frac(loc_y)
    MixFrac[n] = bilger(Y_C[n], Y_H[n], Y_O[n])
    D[n] = gas.thermal_conductivity / (rho[n]*gas.cp_mass)
    
dZdn = df_upwind(MixFrac, x, u)
kai = np.zeros(N)
for n in range(N):
    kai[n] = 2 * D[n] * pow(dZdn[n], 2)

fn = os.path.split(raw_data_path)[-1].split('.txt')[0][:-3] + "transformed.txt"
fp = os.path.join(output_dir, fn)

fout = open(fp, 'w')
header = "{:>18s}".format("rho")
header += "{:>18s}{:>18s}{:>18s}{:>18s}".format("Y_C", "Y_H", "Y_O", "Y_N")
header += "{:>18s}{:>18s}{:>18s}{:>18s}\n".format('Z', 'dZdn', 'D', 'kai')
fout.write(header)
for n in range(N):
    content = '{:>18.8e}'.format(rho[n])
    content += '{:>18.8e}{:>18.8e}{:>18.8e}{:>18.8e}'.format(Y_C[n], Y_H[n], Y_O[n], Y_N[n])
    content += '{:>18.8e}{:>18.8e}{:>18.8e}{:>18.8e}\n'.format(MixFrac[n], dZdn[n], D[n], kai[n])
    fout.write(content)
fout.close()
