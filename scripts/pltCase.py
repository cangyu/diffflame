import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import re

if len(sys.argv) != 2:
    print("Usage: python3 pltCase.py data.txt")
    exit(-1)
else:
    data_file_path = sys.argv[1]
    if not os.path.exists(data_file_path):
        print("Not found!")
        exit(-2)

data = np.loadtxt(data_file_path, skiprows=1)
N, _ = data.shape

mdot_str = re.findall(r'\d+\.?\d*e?[-+]?\d+', data_file_path)
mf = float(mdot_str[0])
mo = float(mdot_str[1])
mL = mf
mR = -mo

z = data[:, 0]
u = data[:, 1]
V = data[:, 2]
T = data[:, 3]
Y_H2 = data[:, 5]
Y_O2 = data[:, 8]
Y_H2O = data[:, 10]
Y_CH4 = data[:, 18]
Y_CO = data[:, 19]
Y_CO2 = data[:, 20]
Y_NO = data[:, 40]
Y_NO2 = data[:, 41]
Y_N2 = data[:, 52]
Y_AR = data[:, 53]


fig1 = plt.figure()
Yax = fig1.add_subplot(1,1,1)
Yax.set_title("mf={} mo={}".format(mf, mo))
Yax.plot(z, 10*Y_AR, label='Y_AR x 10')
Yax.plot(z, Y_CH4, label='Y_CH4')
Yax.plot(z, Y_H2, label='Y_H2')
Yax.plot(z, Y_N2, label='Y_N2')
Yax.plot(z, Y_O2, label='Y_O2')
Yax.plot(z, 10*Y_CO2, label='Y_CO2 x 10')
Yax.plot(z, 10*Y_CO, label='Y_CO x 10')
Yax.plot(z, 1e4*Y_NO2, label='Y_NO2 x 1e4')
Yax.plot(z, 1e3*Y_NO, label='Y_NO x 1e3')
Yax.set_ylabel('Y')
Yax.legend()
Tax = Yax.twinx()
Tax.plot(z, T, 'r')

plt.show()
