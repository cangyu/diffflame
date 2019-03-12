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

data = np.loadtxt(task_file_path)
n = len(data)

def counterflow(x):
    mf = x[0]
    mo = x[1]
    counterflow=subprocess.Popen(["../main.out", "{:f}".format(mf), " {:f}".format(mo)], cwd='./ChemTbl/')
    counterflow.wait()

p = Pool(os.cpu_count())
p.map(counterflow, data)
    