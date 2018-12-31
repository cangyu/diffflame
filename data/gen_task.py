import numpy as np

mf_left, mf_right = map(float, input('[mf_begin, mf_end):').split())
step = float(input('Step: '))
ratio = float(input('O/F ratio: '))

mf = np.arange(mf_left, mf_right, step)

fout = open('task.txt', 'w')
for e in mf:
    fout.write('{:e}\t{:e}\n'.format(e, e*ratio))
fout.close()
