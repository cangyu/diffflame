import os
import sys
import subprocess
from multiprocessing import Pool
import numpy as np
import re

np = os.cpu_count()

data_dir = os.path.join('..', 'data')
if not os.path.exists(data_dir):
    print("Input data directory not found!")
    exit(-2)
output_dir = data_dir

data_list = os.listdir(data_dir)

def transform(f):
    if f.endswith('_raw.txt'):
        case_name = f[:-8]
        print('Processing {} ...'.format(case_name))
        data_file_path = os.path.join(data_dir, f)
        app = subprocess.Popen(["python3", "phy2z.py", data_file_path, output_dir])
        app.wait()

p = Pool(np)
p.map(transform, data_list)
