import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

chemtbl_dir = os.path.join('..', 'ChemTbl')
chemtbl_path = os.path.join(chemtbl_dir, 'ChemTbl.txt')

data = np.loadtxt(chemtbl_path, skiprows=1)
kai = np.log10(data[:, 1])
Z = data[:, 2]*0.01
T = data[:, -2]

fig = plt.figure()
ax = fig.add_subplot(111)
T_bar_norm = colors.Normalize(vmin=300, vmax=2400)
tbl_pic = ax.scatter(kai, Z, c=T, s=10, cmap=plt.cm.rainbow, norm=T_bar_norm, alpha=0.8)
ax.set_xlabel('log10(kai_st)')
ax.set_ylabel('Z')

T_color_bar = plt.colorbar(tbl_pic)
T_color_bar.set_label('T/K')

plt.show()
