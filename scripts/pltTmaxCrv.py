import os
import sys
import numpy as np
import math
from matplotlib import pyplot as plt
import re
import xlwt


def lin_ratio(a, b, x):
    return (x - a) / (b - a)


def relaxation(a, b, alpha):
    return (1 - alpha) * a + alpha * b


if len(sys.argv) != 2:
    print("Usage: python3 pltTmaxCrv.py data_dir")
    exit(-1)
else:
    data_path = sys.argv[1]
    if not os.path.exists(data_path):
        print("Input data directory not found!")
        exit(-2)


Z_st = 0.0496

T_stat = []
kai_stat = []

T_ignition = []
m_f = []
m_o = []
kai_ignition = []


data_set = {}
for f in os.listdir(data_path):
    fn = f.split('.txt')[0]
    if fn.endswith('_raw'):
        case_name = fn[:-4]
    if fn.endswith('_transformed'):
        case_name = fn[:-12]
    if case_name not in data_set:
        data_set[case_name] = 0
    data_set[case_name] += 1

paired_set = []
for case_name in data_set:
    if data_set[case_name] == 2:
        paired_set.append(case_name)

for case_name in paired_set:
    fraw = case_name + '_raw.txt'
    ftrans = case_name + '_transformed.txt'

    # Load transformed data
    fp = os.path.join(data_path, ftrans)
    data = np.loadtxt(fp, skiprows=1)
    kai = data[:, 8]
    Z = data[:, 5]
    # Load raw data
    fp = os.path.join(data_path, fraw)
    data = np.loadtxt(fp, skiprows=1)
    T = data[:, 3]

    # Filter rows where kai is too small
    col_rm = []
    n = len(kai)
    for i in range(n):
        if kai[i] < 1e-20:
            col_rm.append(i)

    kai = np.delete(kai, col_rm)
    Z = np.delete(Z, col_rm)
    T = np.delete(T, col_rm)

    # Find kai where Z = Z_st
    n = len(kai)
    upper_idx = 0
    lower_idx = 0

    err = 10.0
    for i in range(n):
        if Z[i] > Z_st:
            local_err = Z[i] - Z_st
            if local_err < err:
                err = local_err
                upper_idx = i

    err = 10.0
    for i in range(n):
        if Z[i] < Z_st:
            local_err = Z_st - Z[i]
            if local_err < err:
                err = local_err
                lower_idx = i

    r = lin_ratio(Z[lower_idx], Z[upper_idx], Z_st)
    kai_st = relaxation(kai[lower_idx], kai[upper_idx], r)
    kai_stat.append(kai_st)

    # Find max T
    cur_Tmax = max(T)
    T_stat.append(cur_Tmax)
    if cur_Tmax > 350:
        r = re.split(r'[_ \n]', case_name)
        cur_mf = float(r[0][3:])
        cur_mo = float(r[1][3:])
        m_f.append(cur_mf)
        m_o.append(cur_mo)
        T_ignition.append(cur_Tmax)
        kai_ignition.append(kai_st)

    print('{:<24s} T_max={:>8.2f}K kai_st={:>10.2f}'.format(case_name, cur_Tmax, kai_st))

# Save all solution
n = len(kai_stat)
with open('all_solution.txt', 'w') as f:
    for i in range(n):
        f.write('{:e}\t{:e}\n'.format(kai_stat[i], T_stat[i]))

# Save ignition solution
book = xlwt.Workbook(encoding='utf-8', style_compression=0)
sheet = book.add_sheet('solution', cell_overwrite_ok=True)
sheet.write(0, 0, 'm_f')
sheet.write(0, 1, 'm_o')
sheet.write(0, 2, 'T_max')
sheet.write(0, 3, 'kai_st')
n = len(T_ignition)
for k in range(n):
    sheet.write(k + 1, 0, m_f[k])
    sheet.write(k + 1, 1, m_o[k])
    sheet.write(k + 1, 2, T_ignition[k])
    sheet.write(k + 1, 3, kai_ignition[k])
book.save('ignition_solution.xls')

# Plot
plt.scatter(np.log10(kai_stat), T_stat)
plt.xlabel('log10(kai_st)')
plt.ylabel('Tmax')
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig('S-Curve.png', dpi=300)
