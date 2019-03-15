import os
import sys
import math
import numpy as np
from matplotlib import pyplot as plt

base_data = np.loadtxt('../ChemTbl/sCurve_grimech.dat', skiprows=1)

kai_st1 = base_data[:, 0]
Tmax1 = base_data[:, 1]

our_data = np.loadtxt('../ChemTbl/all_solution.txt')
kai_st2 = our_data[:, 0]
Tmax2 = our_data[:, 1]

plt.scatter(np.log10(kai_st1), Tmax1, label='RefSol')
plt.scatter(np.log10(kai_st2), Tmax2, label='PhySol')
plt.xlabel('log10(kai_st)')
plt.ylabel('Tmax/K')
plt.legend()

plt.show()
