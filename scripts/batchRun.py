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
        print('Case for mf={} mo={} already exists!'.format(mf, mo))
        dup_exist = True

    if (not dup_exist) or (dup_exist and OverwriteExisting):
        counterflow=subprocess.Popen(["../src/main.out", "{}".format(mf), "{}".format(mo), "{}".format(domain_length)], cwd=output_dir)
        counterflow.wait()

p = Pool(NumOfProc)
p.map(counterflow, data)
