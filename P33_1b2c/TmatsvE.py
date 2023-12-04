#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math

plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14,6))

# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

# read t-matrix data from file
data = np.loadtxt('dataInf_orig.in', skiprows=1)
E_data = data[:,0]
Treal_data = data[:,5]
dTreal_data = data[:,6]
Timag_data = data[:,7]
dTimag_data = data[:,8]


# output = np.loadtxt(f"data/tmatrix_fit{fitNum}.out", skiprows=1)



# ---------------------------Plot tmatrix---------------------------
axs[0].errorbar(E_data, Treal_data, yerr=dTreal_data, fmt='.',
              ecolor='red', color='red', capsize=2, zorder=5.0)

axs[0].set_ylabel('Re $T$')


axs[1].errorbar(E_data, Timag_data, yerr=dTimag_data, fmt='.',
              ecolor='blue', color='blue', capsize=2, zorder=5.0)

axs[1].set_ylabel('Im $T$')


for ax in axs:
    ax.set_xlim(1.08, max(E_data)*1.05)
    ax.set_ylim(-0.6, 1.0)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # xtickloc = []
    # ax.xticks(xtickloc)
    ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
    ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
    ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
    ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')



plt.show()
