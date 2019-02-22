import os
import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data_dir = os.path.join('.', 'ChemTbl')
data = []

for f in os.listdir(data_dir):
    if f.endswith('txt'):
        fp = os.path.join(data_dir, f)
        cur_data = list(np.loadtxt(fp, skiprows=1))
        data = data + cur_data

n = len(data)
kzt = []
for i, e in enumerate(data):
    kai = e[1]
    Z = e[2]
    T = e[-2]
    if kai >= 1e-10:
        kzt.append([math.log10(kai), Z, T])

kzt = np.array(kzt)

# with open('ChemTbl.txt','w') as f:
#     np.savetxt(f, kzt, fmt='%e')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('log10(Kai)')
ax.set_ylabel('Mixture Fraction')
ax.set_zlabel('Temperature')
ax.scatter(kzt[:, 0], kzt[:, 1], kzt[:, 2])

plt.show()
