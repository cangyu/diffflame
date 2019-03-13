import os
import sys
import subprocess
from multiprocessing import Pool
import numpy as np

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
    cur_case_output_name = 'mf={}_mo={}_raw.txt'.format(mf, mo)
    if not os.path.exists(os.path.join(output_dir, cur_case_output_name)):
        counterflow=subprocess.Popen(["../src/main.out", "{:f}".format(mf), " {:f}".format(mo)], cwd=output_dir)
        counterflow.wait()
    else:
        print('Case for mf={} mo={} already exists!'.format(mf, mo))

p = Pool(os.cpu_count()//2)
p.map(counterflow, data)
