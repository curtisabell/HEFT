#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

sys.path.append('../heftCode')
import readHEFTConfig
HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()
n_ch = HEFT.n_ch
n_bare = HEFT.n_bare

# -----------------------Stuff for saving figs----------------------
saveFigs = False
figType = 'png'

cmdArgs = sys.argv
if (len(cmdArgs)>1):
    saveArg = str(cmdArgs[1])
    if (saveArg.lower()=='save'):
        saveFigs = True
if (len(cmdArgs)>2):
    typeArg = str(cmdArgs[2])
    if (typeArg.lower()=='pdf'):
        figType = 'pdf'


plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, ax = plt.subplots()



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


output = np.loadtxt(f"data/tmatrix_fit{fitNum}.out")
E_output = output[:,0]
Treal_output = output[:,1]
Timag_output = output[:,2]


# ---------------------------Plot tmatrix---------------------------
plot_real, = ax.plot(E_output, Treal_output, color='red', linestyle='-'
                    , linewidth='1.7', label='Treal')
ax.errorbar(E_data, Treal_data, yerr=dTreal_data, fmt='.',
              ecolor='red', color='red', capsize=2, zorder=5.0)


plot_imag, = ax.plot(E_output, Timag_output, color='blue', linestyle='-'
                    , linewidth='1.7', label='Timag')
ax.errorbar(E_data, Timag_data, yerr=dTimag_data, fmt='.',
              ecolor='blue', color='blue', capsize=2, zorder=5.0)


ax.plot([E_data[0],E_data[-1]], [0,0], color='black', linewidth=0.5)

ax.legend([plot_real, plot_imag], ['Re $T$', 'Im $T$'], prop={'size': 18}
          , loc='upper right', frameon=False, handletextpad=0.1)
# ax.set_title(f'P33 1b1c, fit {fitNum}')
ax.set_xlabel('E (GeV)')
ax.set_xlim(min(E_output), max(E_output))
ax.set_ylim(-0.6, 1.05)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
# xtickloc = []
# ax.xticks(xtickloc)
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')


outfile = f'figs/{int(n_bare)}b{int(n_ch)}c_Tmat_fit{fitNum}.{figType}'
if saveFigs:
    savefig(outfile, bbox_inches='tight')
else:
    plt.show()
