#!/usr/bin/python3
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time

nPoints = [64, 128, 256, 512]
nCores = [1, 2, 4, 8, 16]
t = np.zeros((len(nCores), len(nPoints)))

calcSpeedup = False

if calcSpeedup:
    with open('makefile', 'r') as f:
        makefile_lines = f.readlines()
    with open('HEFTFinite.config', 'r') as f:
        HEFT_lines = f.readlines()

    for iPoint, points in enumerate(nPoints):
        print(f'nPoints = {points}')

        # Change the number of points to run with
        new_points_line = f'L_points   {points}\n'
        HEFT_lines[8] = new_points_line
        with open('HEFTFinite.config', 'w') as f:
            f.writelines(HEFT_lines)

        print('   nCores = ', end='', flush=True)
        for iCore, core in enumerate(nCores):
            # Print progress
            print(f'{core}', end='', flush=True)
            if iCore == len(nCores)-1:
                print('')
            else:
                print(', ', end='', flush=True)

            # Change the number of cores to compile with
            new_coarr_line = f'CFLAGS   = -O2 -mkl -real-size 64 -coarray-num-images={core}\n'
            makefile_lines[3] = new_coarr_line
            with open('makefile', 'w') as f:
                f.writelines(makefile_lines)

            # Run and time the finite volume HEFT code
            sp.run('make clean', shell=True, stdout=sp.DEVNULL)
            sp.run('make fin.x', shell=True, stdout=sp.DEVNULL)
            t1 = time.perf_counter()
            sp.run('./fin.x', shell=True, stdout=sp.DEVNULL)
            t2 = time.perf_counter()
            t[iCore, iPoint] = t2 - t1
    print('')

    speedup = np.zeros_like(t)
    for i in range(len(nPoints)):
        speedup[:, i] = max(t[:, i]) / t[:, i]
    np.savetxt('speedup.data', speedup)
else:
    speedup = np.loadtxt('speedup.data')

# ------------------------------------------------------------------

df = pd.DataFrame(data=speedup,
                  index=[f'{iCore}' for iCore in nCores],
                  columns=[f'N={iPoint}' for iPoint in nPoints])

df.plot()
ax = plt.gca()
ax.set(xlabel='nCores', ylabel='Speedup')
# df.plot(x='nCores', y='Speedup')
plt.show()
