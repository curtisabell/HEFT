#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

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

# ------------------------------------------------------------------
plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, ax = plt.subplots()

# ------------------------------------------------------------------
corr_filename = f'data/correlation_lqcd_fit{fitNum}.out'
data_corr = np.loadtxt(corr_filename)
nMass = int(len(data_corr[:,0]) / 2)
nt = int(len(data_corr[0,:])) - 1

m_pi2 = np.zeros(nMass)
L = np.zeros(nMass)
t = np.zeros((nMass,nt))
Gt = np.zeros((nMass,nt))

for i in range(nMass):
    m_index = i*2
    L_index = i*2 + 1

    m_pi2[i] = data_corr[m_index,0]
    L[i] = data_corr[L_index,0]
    t[i,:] = data_corr[m_index,1:]
    Gt[i,:] = data_corr[L_index,1:]

colours = ['red', 'blue']

# choose which m_pi2 to plot
m_pi2_choice = 0.265
ind = np.argmin(np.abs(m_pi2 - m_pi2_choice))

ax.plot(t[ind,:], Gt[ind,:], color='red')


ax.set_title(f'S11 {int(n_bare)}b{int(n_ch)}c, L = {L[ind]:.1f} fm, ' + '$m_{\pi}^{2}$' + f' = {m_pi2[ind]:.3f} GeV' + '${}^2$')
ax.set_xlabel('$t$ (fm)')
ax.set_ylabel('$G(t)$')

outfile = f'figs/S11_{int(n_bare)}b{int(n_ch)}c_correlation_L{int(L[ind]*10)}_mpisq{int(m_pi2[ind]*1000)}_fit{fitNum}.{figType}'


# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

if saveFigs:
    savefig(outfile, bbox_inches='tight')
else:
    plt.show()
