import os
import numpy as np
import math
from matplotlib import pyplot as plt
import re
import xlwt

data_path = os.path.join('.', 'ChemTbl')


def lin_ratio(a, b, x):
    return (x - a) / (b - a)


def relaxation(a, b, alpha):
    return (1 - alpha) * a + alpha * b


Z_st = 0.0496

T_stat = []
kai_stat = []

T_ignition = []
m_f = []
m_o = []
kai_ignition = []

'''
cur_task = list(np.loadtxt('task.txt'))
n1 = len(cur_task)
print(n1)

todo_pair = []
for e in cur_task:
    todo_pair.append((e[0], e[1]))

solved_pair = []
for f in os.listdir(data_path):
    r = re.split(r'[_ \n]', f[:-4])
    cur_mf = float(r[0][3:])
    cur_mo = float(r[1][3:])
    solved_pair.append((cur_mf, cur_mo))
n2 = len(solved_pair)
print(n2)

filtered_pair = list(set(todo_pair).difference(set(solved_pair)))
n3 = len(filtered_pair)
print(n3)
with open('filtered_task.txt', 'w') as f:
    for e in filtered_pair:
        f.write('{:e}\t{:e}\n'.format(e[0], e[1]))
'''

for f in os.listdir(data_path):
    # Load original data
    fp = os.path.join(data_path, f)
    data = np.loadtxt(fp, skiprows=1, usecols=[1, 2, 9])
    kai = data[:, 0]
    Z = data[:, 1]
    T = data[:, 2]

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
        r = re.split(r'[_ \n]', f[:-4])
        cur_mf = float(r[0][3:])
        cur_mo = float(r[1][3:])
        m_f.append(cur_mf)
        m_o.append(cur_mo)
        T_ignition.append(cur_Tmax)
        kai_ignition.append(kai_st)

    print('{:<36s} T_max={:>8.2f}K kai_st={:>10.2f}'.format(f[:-4], cur_Tmax, kai_st))

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
