import os
import math
import numpy as np


def linear_interp(a, b, t: float):
    return (1-t) * a + t * b


def convert_to_list(sol):
    I, J, _ = sol.shape
    rl = []
    for i in range(I):
        for j in range(J):
            rl.append(sol[i][j])
    return rl


NumOfChemTblRcrdElem = 11
KAI_IDX = 1
Z_IDX = 2
CH4_IDX = -7
CO_IDX = -6
H2_IDX = -5
H2O_IDX = -4
CO2_IDX = -3
T_IDX = -2
NO_IDX = -1


def transform_z(e):
    e[Z_IDX] = round(100 * e[Z_IDX])


def filter_kai_and_trans_z(rl):
    fdata = []
    for e in rl:
        kai = e[KAI_IDX]
        if not math.isclose(kai, 0.0, abs_tol=1e-10):
            transform_z(e)
            fdata.append(e)
    return np.array(fdata)


def interp_z(rl):
    interp_record = []
    n = len(rl)
    for i in range(1, n):
        left_z = int(rl[i-1][Z_IDX])
        right_z = int(rl[i][Z_IDX])
        gap = abs(right_z - left_z)
        for j in range(1, gap):
            cur_interp_rec = linear_interp(rl[i-1], rl[i], j/gap)
            interp_record.append(cur_interp_rec)
    return list(rl) + interp_record


def filter_z(rl):
    rl.sort(key=lambda x: x[Z_IDX])
    n = len(rl)

    rg = []
    cur_z = -1.0
    for i in range(n):
        if rl[i][Z_IDX] != cur_z:
            rg.append(i)
            cur_z = rl[i][Z_IDX]
    rg.append(n)

    ret = []
    for i in range(1, len(rg)):
        l = rg[i-1]
        r = rg[i]
        e = sum(rl[l:r])/(r-l)
        ret.append(e)

    return np.array(ret)


def add_missing_z(rl):
    ret = []

    z_head = int(rl[0][Z_IDX])
    if z_head > 0:
        for i in range(z_head):
            t = np.copy(rl[0])
            t[Z_IDX] = i*1.0
            ret.append(t)

    z_tail = int(rl[-1][Z_IDX])
    if z_tail < 100:
        for i in range(z_tail+1, 101):
            t = np.copy(rl[-1])
            t[Z_IDX] = i*1.0
            ret.append(t)

    full_list = ret + list(rl)
    full_list.sort(key=lambda x: x[Z_IDX])

    return np.array(full_list)


def add_extinction_and_do_average(rl):
    k_num = 75  # 按区间统计
    z_num = 100  # 按点统计

    lkl = -10.0
    lkr = 4.0
    lkg = lkr - lkl

    def calc_k_idx(lk):
        return int((lk - lkl)/lkg * k_num)

    def calc_kai(i):
        return math.pow(10, linear_interp(lkl, lkr, (i+0.5)/k_num))

    sol = np.zeros((k_num, z_num+1, NumOfChemTblRcrdElem))
    cold_sol = [[[] for j in range(z_num+1)] for i in range(k_num)]
    hot_sol = [[[] for j in range(z_num+1)] for i in range(k_num)]

    for k, e in enumerate(rl):
        lk = math.log10(e[KAI_IDX])
        z = e[Z_IDX]
        T = e[T_IDX]
        if lkl < lk < lkr:
            i = calc_k_idx(lk)
            j = int(z)
            if T < 310:
                cold_sol[i][j].append(k)
            else:
                hot_sol[i][j].append(k)

    for i in range(k_num):
        kai = calc_kai(i)
        for j in range(z_num+1):
            if (not hot_sol[i][j]) and (not cold_sol[i][j]):
                z_u = 0.01*j
                sol[i][j] = np.array([1e-4, kai, j, 0.0, 8/9*z_u, 0.0, 1/9*z_u, 0.0, 0.0, 300.0, 0.0]) #燃料侧是甲烷-氢气1:1
            elif (not hot_sol[i][j]) and cold_sol[i][j]:
                sol[i][j] = sum([rl[k] for k in cold_sol[i][j]])/len(cold_sol[i][j])
            else:
                sol[i][j] = sum([rl[k] for k in hot_sol[i][j]])/len(hot_sol[i][j])
            
            sol[i][j][KAI_IDX] = kai

    for j in range(z_num+1):
        hot_stat = []
        for i in range(k_num):
            if sol[i][j][T_IDX] > 310:
                hot_stat.append(i)
        
        if hot_stat:
            y_ch4 = min([sol[i][j][CH4_IDX] for i in hot_stat])
            y_co = min([sol[i][j][CO_IDX] for i in hot_stat])
            y_h2 = min([sol[i][j][H2_IDX] for i in hot_stat])
            y_h2o = max([sol[i][j][H2O_IDX] for i in hot_stat])
            y_co2 = max([sol[i][j][CO2_IDX] for i in hot_stat])
            T = max([sol[i][j][T_IDX] for i in hot_stat])
            y_no = min([sol[i][j][NO_IDX] for i in hot_stat])

            for i in range(k_num):
                kai = calc_kai(i)
                if kai < 1.0 and sol[i][j][T_IDX] <= 310:
                    sol[i][j][CH4_IDX] = y_ch4
                    sol[i][j][CO_IDX] = y_co
                    sol[i][j][H2_IDX] = y_h2
                    sol[i][j][H2O_IDX] = y_h2o
                    sol[i][j][CO2_IDX] = y_co2
                    sol[i][j][T_IDX] = T
                    sol[i][j][NO_IDX] = y_no

    for i in range(k_num):
        kai = calc_kai(i)
        sol[i][0] = np.array([1e-4, kai, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 300.0, 0.0])

    return sol


def output1(rl, fn):
    with open(fn, 'w') as f:
        for tag in ['D', 'kai', 'Z', 'Y_CH4', 'Y_CO', 'Y_H2', 'Y_H2O', 'Y_CO2', 'T', 'Y_NO']:
            f.write('{:>18s}\t'.format(tag))
        f.write('\n')

        for e in rl:
            for i in [0, 1]:
                f.write('{:>18.8e}\t'.format(e[i]))
            mix_frac = int(e[Z_IDX])
            f.write('{:>18d}\t'.format(mix_frac))
            for i in [4, 5, 6, 7, 8, 9, 10]:
                f.write('{:>18.8e}\t'.format(e[i] if e[i] > 0 else 0.0))
            f.write('\n')


def output2(rl, fn):
    I, J, _ = rl.shape
    f = open(fn, 'w')
    for _ in range(15):
        for i in range(I):
            for spec in range(7):
                for j in range(J):
                    kai = rl[i][j][KAI_IDX]
                    Z = int(rl[i][j][Z_IDX])
                    val = max(rl[i][j][spec+4], 0.0)
                    dts = '{:>18.8e}\t{:>18.8e}\t{:>18d}\t{:>18d}\t{:>18.8e}\n'.format(1e-4, kai, spec+1, Z, val)
                    f.write(dts)
    f.close()


if __name__ == '__main__':
    origin_chem_tbl_path = os.path.join('.', 'ChemTbl')
    interp_chem_tbl_path = os.path.join('.', 'Interpolated')
    
    if not os.path.exists(interp_chem_tbl_path):
        os.makedirs(interp_chem_tbl_path)

    data = []
    for f in os.listdir(origin_chem_tbl_path):
        if f.endswith('.txt'):
            d0 = np.loadtxt(os.path.join(origin_chem_tbl_path, f), skiprows=1)
            d1 = filter_kai_and_trans_z(d0)
            d2 = interp_z(d1)
            d3 = filter_z(d2)
            d4 = add_missing_z(d3)
            output1(d4, os.path.join(interp_chem_tbl_path, f))
            data = data + list(d4)

    d5 = add_extinction_and_do_average(data)
    output2(d5, os.path.join(interp_chem_tbl_path, 'chazi3DNOnnnoo.dat'))
    d6 = convert_to_list(d5)
    output1(d6, os.path.join(interp_chem_tbl_path, 'ChemTbl.txt'))
