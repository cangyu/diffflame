import sys
import math
import numpy as np
import cantera as ct
from matplotlib import pyplot as plt


rel_tol = 1e-5
abs_tol = math.sqrt(sys.float_info.epsilon)


def perturbation_delta(x: float):
    return rel_tol * x + abs_tol


def relaxation(a, b, x: float):
    return (1 - x) * a + x * b


def df_central(f, x):
    pnt_num = len(x)
    ret = np.zeros(pnt_num)
    ret[0] = ((x[1] - x[0]) / (x[2] - x[0]) * (f[2] - f[0]) - (x[2] - x[0]) / (x[1] - x[0]) * (f[1] - f[0])) / (x[1] - x[2])
    for i in range(1, pnt_num - 1):
        ret[i] = ((x[i] - x[i - 1]) / (x[i + 1] - x[i]) * (f[i + 1] - f[i]) - (x[i + 1] - x[i]) / (x[i] - x[i - 1]) * (f[i - 1] - f[i])) / (x[i + 1] - x[i - 1])
    ret[-1] = ((x[-1] - x[-3]) / (x[-1] - x[-2]) * (f[-2] - f[-1]) - (x[-1] - x[-2]) / (x[-1] - x[-3]) * (f[-3] - f[-1])) / (x[-3] - x[-2])
    return ret


def df_upwind(f, x, u):
    pnt_num = len(x)
    ret = np.zeros(pnt_num)
    ret[0] = (f[1] - f[0]) / (x[1] - x[0])
    for i in range(1, pnt_num - 1):
        s = -1 if u[i] > 0 else 1
        ret[i] = (f[i + s] - f[i]) / (x[i + s] - x[i])
    ret[-1] = (f[-1] - f[-2]) / (x[-1] - x[-2])
    return ret


def ddf(f, x):
    pnt_num = len(x)
    ret = np.zeros(pnt_num)
    ret[0] = 2.0 / (x[2] - x[1]) * ((f[2] - f[0]) / (x[2] - x[0]) - (f[1] - f[0]) / (x[1] - x[0]))
    for i in range(1, pnt_num - 1):
        ret[i] = 2.0 / (x[i + 1] - x[i - 1]) * ((f[i - 1] - f[i]) / (x[i] - x[i - 1]) + (f[i + 1] - f[i]) / (x[i + 1] - x[i]))
    ret[-1] = 2.0 / (x[-3] - x[-2]) * ((f[-2] - f[-1]) / (x[-1] - x[-2]) - (f[-3] - f[-1]) / (x[-1] - x[-3]))
    return ret


gas = ct.Solution('gri30.xml')
Le = 1.0
gas.transport_model = 'UnityLewis'
K = gas.n_species  # Num of species
NAME = gas.species_names
MW = gas.molecular_weights  # Kg/Kmol

ich4 = gas.species_index('CH4')
ih2 = gas.species_index('H2')
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

P = ct.one_atm # Pa
mf = 0.1  # kg / (m^2 * s)
mo = 0.3  # kg / (m^2 * s)
L = 0.05  # m


def gen_data_str(mf, mo, L):
    c = "mf={}_mo={}_L={}".format(mf, mo, L)
    raw_str = c + "_raw.txt"
    trans_str = c + "_transformed.txt"
    return raw_str, trans_str

raw_data_str, trans_data_str = gen_data_str(mf, mo, L)

raw_data = np.loadtxt(raw_data_str, skiprows=1)
trans_data = np.loadtxt(trans_data_str, skiprows=1)
N = len(raw_data)  # Num of points
assert(len(trans_data)==N), "Transfomed data not consistent with raw data!"

z = raw_data[:, 0]
mL = mf
mR = -mo

CUR = 0
NEXT = 1

rho = np.zeros([2, N])  # kg/m^3
u = np.zeros([2, N])  # m / s
V = np.zeros([2, N])  # s^-1
A = np.zeros([2, N])  # The eigenvalue, kg / (m^3 * s^2)
T = np.zeros([2, N])  # K
Y = np.zeros([2, K, N])  # Mass fraction

mu = np.zeros([2, N])  # Viscosity, Pa * s = Kg / (m * s)
Cp = np.zeros([2, N])  # Specific heat, J / (Kg * K)
Lambda = np.zeros([2, N])  # Thermal conductivity, W / (m * K)
D = np.zeros([2, N])  # Diffusion coefficients, m ^ 2 / s

RR = np.zeros([2, K, N])  # Chemical reaction rate, Kg / (m ^ 3 * s)
RS = np.zeros([2, N])  # Energy source due to chemical reaction, J / (m ^ 3 * s)

dVdz = np.zeros([2, N])
dTdz = np.zeros([2, N])
dYdz = np.zeros([2, K, N])
ddVddz = np.zeros([2, N])
ddTddz = np.zeros([2, N])
ddYddz = np.zeros([2, K, N])

C = K + 4 # Num of unknown per node
Q = C * N # Total num of unknown
phi = np.zeros([2, Q]) # Solution vector
F = np.zeros([2, Q]) # Residual vector
J = np.zeros([Q, Q]) # Jacobian matrix


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


def calc_residual(level):
    # Update properties
    for i in range(N):
        gas.TPY = T[level, i], P, Y[level, :, i]
        mu[level, i] = gas.viscosity
        Lambda[level, i] = gas.thermal_conductivity
        Cp[level, i] = gas.cp_mass
        D[level, i] = Lambda[level, i] / (rho[level, i] * Cp[level, i] * Le)
        w = gas.net_production_rates  # kmol / (m ^ 3 * s)
        RR[level, :, i] = w * MW  # Kg / (m ^ 3 * s)
        h = gas.standard_enthalpies_RT * T[level, i] * ct.gas_constant  # J / Kmol
        RS[level, i] = np.dot(w, h)  # J / (m ^ 3 * s)

    # Compute derivatives
    dVdz[level, :] = df_upwind(V[level, :], z, u[level, :])
    dTdz[level, :] = df_upwind(T[level, :], z, u[level, :])
    for k in range(K):
        dYdz[level, k, :] = df_upwind(Y[level, k, :], z, u[level, :])
    ddVddz[level, :] = ddf(V[level, :], z)
    ddTddz[level, :] = ddf(T[level, :], z)
    for k in range(K):
        ddYddz[level, k, :] = ddf(Y[level, k, :], z)

    # Calculate residuals
    cnt = 0
    for i in range(N):
        F[level, cnt] = rho[level, i] * Cp[level, i] * u[level, i] * dTdz[level, i] - Lambda[level, i] * ddTddz[level, i] + RS[level, i]
        cnt += 1
        for k in range(K):
            F[level, cnt] = rho[level, i] * u[level, i] * dYdz[level, k, i] - D[level, k, i] * ddYddz[level, k, i] - RR[level, k, i]
            cnt += 1
        F[level, cnt] = rho[level, i] * u[level, i] * dVdz[level, i] + rho[level, i] * V[level, i] ** 2 + A[level, i] - mu[level, i] * ddVddz[level, i]
        cnt += 1


'''Initialization'''
for p in range(N):
    u[CUR][p] = raw_data[p, 1]
    V[CUR][p] = raw_data[p, 2]
    T[CUR][p] = raw_data[p, 3]
    A[CUR][p] = raw_data[p, 4]
    for k in range(K):
        Y[CUR][k][p] = raw_data[p, k+5]
    rho[CUR][p] = trans_data[p, 0]
    D[CUR, p] = trans_data[p, 7]

dVdz[CUR, :] = df_upwind(V[CUR, :], z, u[CUR, :])
ddVddz[CUR, :] = ddf(V[CUR, :], z)
dTdz[CUR, :] = df_upwind(T[CUR, :], z, u[CUR, :])
ddTddz = ddf(T[CUR, :], z)
for k in range(K):
    dYdz[CUR, k, :] = df_upwind(Y[CUR, k, :], z, u[CUR, :])
    ddYddz[CUR, k, :] = ddf(Y[CUR, k, :], z)


'''Damped Newton Method'''
calc_residual(CUR)
for j in range(Q):
    rho[NEXT, :] = rho[CUR, :]
    u[NEXT, :] = u[CUR, :]
    V[NEXT, :] = V[CUR, :]
    A[NEXT, :] = A[CUR, :]
    T[NEXT, :] = T[CUR, :]
    Y[NEXT, :, :] = Y[CUR, :, :]
    # Add perturbation
    node_idx = int(j / (C - 1))
    var_idx = j - node_idx * C
    if var_idx is 0:
        delta = perturbation_delta(rel_tol, abs_tol, T[CUR, node_idx])
        T[NEXT, node_idx] += delta
    elif var_idx is C - 1:
        delta = perturbation_delta(rel_tol, abs_tol, A[CUR, node_idx])
        A[NEXT, node_idx] += delta
    else:
        spec_idx = var_idx - 1
        delta = perturbation_delta(rel_tol, abs_tol, Y[CUR, spec_idx, node_idx])
        Y[NEXT, spec_idx, node_idx] += delta

    calc_residual(NEXT)
    J[:, j] = (F[NEXT, :] - F[CUR, :]) / delta

print(J.shape)
