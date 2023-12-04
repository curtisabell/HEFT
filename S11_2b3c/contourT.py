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
m_pi0 = 0.1385

plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, ax = plt.subplots()

# ---------------------------Get fit info---------------------------
# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

# ------------------------Check to save figs------------------------
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


# ------------------------------------------------------------------
data = np.loadtxt(f'data/TMatGrid_fit{fitNum}.out')
E_real = data[0,1:]
E_imag = data[1:,0]
TGrid = data[1:,1:]

cp = ax.contourf(E_real, E_imag, np.log(TGrid))
# fig.colorbar(cp)



# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')



if saveFigs:
    outFile = f'figs/S11_{n_bare}b{n_ch}c_TMatContour_' \
        + f'fit{fitNum}.{figType}'
    savefig(outFile, bbox_inches='tight')

else:
    plt.show()
