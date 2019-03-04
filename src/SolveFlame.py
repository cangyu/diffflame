import sys
import math
import numpy as np
import cantera as ct
from matplotlib import pyplot as plt


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


eps = sys.float_info.epsilon
rel_tol = math.sqrt(eps)
abs_tol = math.sqrt(eps)

gas = ct.Solution('gri30.cti')
K = gas.n_species  # Num of species
NAME = gas.species_names
MW = gas.molecular_weights  # Kg/Kmol

ich4 = gas.species_index('CH4')
ih2 = gas.species_index('H2')
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

P = ct.one_atm
Tin = 300.0  # K
mf = 0.018  # kg / (m^2 * s)
mo = 0.054  # kg / (m^2 * s)
L = 0.1  # m
N = 251  # Num of points
Le = 1.0  # Lewis Number

z = np.linspace(0, L, N)
mL = mf
mR = -mo

CUR = 1
NEXT = 2

rho = np.zeros([2, N])  # Density, kg/m^3
u = np.zeros([2, N])  # m / s
V = np.zeros([2, N])  # s^-1
A = np.zeros([2, N])  # The eigenvalue, kg / (m^3 * s^2)
T = np.zeros([2, N])  # K
Y = np.zeros([2, K, N])  # Mass fraction

mu = np.zeros([2, N])  # Viscosity, Pa * s = Kg / (m * s)
Cp = np.zeros([2, N])  # Specific heat, J / (Kg * K)
Lambda = np.zeros([2, N])  # Thermal conductivity, W / (m * K)
D = np.zeros([2, K, N])  # Diffusion coefficients, m ^ 2 / s

RR = np.zeros([2, K, N])  # Chemical reaction rate, Kg / (m ^ 3 * s)
RS = np.zeros([2, N])  # Energy source due to chemical reaction, J / (m ^ 3 * s)

dVdz = np.zeros([2, N])
dTdz = np.zeros([2, N])
dYdz = np.zeros([2, K, N])
ddVddz = np.zeros([2, N])
ddTddz = np.zeros([2, N])
ddYddz = np.zeros([2, K, N])

C = 1 + K + 1
unknown_num = C * N
phi = np.zeros([2, unknown_num])
F = np.zeros([2, unknown_num])
J = np.zeros([unknown_num, unknown_num])


def build_phi(temp, mass_fraction, eigenvalue, pnt_num, dst):
    cnt = 0
    for i in range(pnt_num):
        dst[cnt] = temp[i]
        cnt = cnt + 1
        for k in range(K):
            dst[cnt] = mass_fraction[k][i]
            cnt = cnt + 1
        dst[cnt] = eigenvalue[i]
        cnt = cnt + 1


def calc_residual(src):
    pass


if __name__ == '__main__':
    '''Simple initialization'''
    for p in range(N):
        cur_ratio = p / (N - 1)
        T[CUR][p] = Tin + (2200.0 - Tin) * math.exp(-(p - N / 2) ** 2 / (2 * 500))
        Y[CUR][ich4][p] = relaxation(8.88370331e-01, 0.0, cur_ratio)
        Y[CUR][ih2][p] = relaxation(1.11629669e-01, 0.0, cur_ratio)
        Y[CUR][io2][p] = relaxation(0.0, 2.33009709e-01, cur_ratio)
        Y[CUR][in2][p] = relaxation(0.0, 7.66990291e-01, cur_ratio)
        gas.TPY = T[CUR][p], P, [Y[CUR][k][p] for k in range(K)]
        rho[CUR][p] = gas.density_mass

    u[CUR][0] = mL / rho[CUR][0]
    u[CUR][-1] = mR / rho[CUR][-1]
    for p in range(1, N - 1):
        cur_ratio = p / (N - 1)
        u[CUR][p] = relaxation(u[CUR][0], u[CUR][-1], cur_ratio)

    V[CUR, :] = -df_central(rho[CUR, :] * u[CUR, :], z) / (2 * rho[CUR, :])
    V[CUR][0] = 0.0
    V[CUR][-1] = 0.0

    dVdz[CUR, :] = df_upwind(V[CUR, :], z, u[CUR, :])
    ddVddz[CUR, :] = ddf(V[CUR, :], z)
    lhs1 = np.dot(rho[CUR, :] * u[CUR, :], dVdz[CUR, :])
    lhs2 = np.dot(rho[CUR, :] * V[CUR, :], V[CUR, :])
    rhs2 = np.dot(mu[CUR, :], ddVddz[CUR, :])
    A[CUR, :] = (rhs2 - lhs1 - lhs2) / N

    '''Time-stepping'''
	#TODO
