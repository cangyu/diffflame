import os
import sys
import subprocess
from multiprocessing import Pool
import numpy as np

NumOfProc = os.cpu_count()
OverwriteExisting = True

if len(sys.argv) != 2:
    print("Usage: python3 batchRun.py task.txt")
    exit(-1)
else:
    task_file_path = sys.argv[1]
    if not os.path.exists(task_file_path):
        print("Not found!")
        exit(-2)

output_dir = os.path.join('..', 'data')
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

data = np.loadtxt(task_file_path)
n = len(data)

def counterflow(x):
    mf = x[0]
    mo = x[1]
    domain_length = x[2]
    cur_case_output_name = 'mf={}_mo={}_L={}_raw.txt'.format(mf, mo, domain_length)

    dup_exist = False
    if os.path.exists(os.path.join(output_dir, cur_case_output_name)):
        dup_exist = True
        info = 'Data of mf={} mo={} L={} already exists!'.format(mf, mo, domain_length)
        if OverwriteExisting:
            info += ' Will be Overwritten!'
        else:
            info += ' Skip!'
        print(info)

    if (not dup_exist) or (dup_exist and OverwriteExisting):
        prog = "../src/main.out"
        arg1 = "{}".format(mf)
        arg2 = "{}".format(mo)
        arg3 = "{}".format(domain_length)
        cmd = [prog, arg1, arg2, arg3]
        counterflow=subprocess.Popen(cmd, cwd=output_dir)
        counterflow.wait()

p = Pool(NumOfProc)
p.map(counterflow, data)
