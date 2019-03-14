import argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('-a', '--append', help='Append mass flux pairs to existing task file')
args = ap.parse_args() 
if args.append:
    fn = args.append
    fout = open(fn, 'a')
else:
    fout = open('task.txt', 'w')

mf_left, mf_right = map(float, input('[mf_begin, mf_end):').split())
step = float(input('Step: '))
ratio = float(input('O/F ratio: '))
L = float(input('Domain length: '))

mf = np.arange(mf_left, mf_right, step)
for e in mf:
    fout.write('{:>18.6e}{:>18.6e}{:>18.6e}\n'.format(e, e*ratio, L))
fout.close()
