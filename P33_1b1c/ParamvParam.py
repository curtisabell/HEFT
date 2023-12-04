#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, ax = plt.subplots()


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


params = np.loadtxt('renormFits_8GeVboth.params', skiprows=2)
m_bare = params[:,0]
gBare = params[:,1]
Lambdas = params[:,2]
vCh = params[:,3]


plotPot = True
if plotPot:
    ax.plot(Lambdas, gBare, 'red', linewidth=1.5, label='$g_{\pi N}^{\Delta}$')
    ax.plot(Lambdas, vCh*1, 'blue', linewidth=1.5, label='$v_{\pi N\pi N}$') #  \\times 5
    ax.plot(Lambdas, Lambdas*0, 'black', linewidth=0.5, label='_nolegend')
    ax.set_ylabel('    ')
else:
    ax.plot(Lambdas, m_bare, 'red', linewidth=1.5, label='$m_{\Delta^{0}}$')
    # ax.set_ylabel('$m_{\Delta^{0}}$ (GeV)')
    ax.set_ylabel('GeV')

plt.legend(loc='best', numpoints=1
           , prop={'size': 18}, frameon=False)

ax.set_xlim([0.6, 8.0])
ax.set_xlabel('$\Lambda$ (GeV)')



# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

if saveFigs:
    outFile = f'figs/1b1c_params.{figType}'
    savefig(outFile, bbox_inches='tight')
else:
    plt.show()
