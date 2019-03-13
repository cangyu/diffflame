import os
import sys
import subprocess
from multiprocessing import Pool
import numpy as np
import re

data_dir = os.path.join('..', 'data')
if not os.path.exists(data_dir):
    print("Input data directory not found!")
    exit(-2)
output_dir = data_dir

data_set = {}
for f in os.listdir(data_dir):
    fn = f.split('.txt')[0]
    
    if fn.endswith('_raw'):
        case_name = fn[:-4]
    if fn.endswith('_transformed'):
        case_name = fn[:-12]
    
    if case_name not in data_set:
        data_set[case_name] = 0
    
    data_set[case_name] += 1


target_file_list = []
for f in os.listdir(data_dir):
    if f.endswith('_raw.txt'):
        case_name = f[:-8]
        if data_set[case_name] == 1:
            data_file_path = os.path.join(data_dir, f)
            target_file_list.append(data_file_path)

for f in target_file_list:
    print('Processing {} ...'.format(os.path.split(f)[-1]))
    transform = subprocess.Popen(["python3", "phy2z.py", f, output_dir])
    transform.wait()
