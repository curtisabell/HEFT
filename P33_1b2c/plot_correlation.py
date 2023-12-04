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
corr_filename = f'data/correlation_fit{fitNum}.out'
with open(corr_filename, 'r') as f:
    line = f.readline().strip().split()
    L = float(line[1])
    line = f.readline().strip().split()
    m_pi = float(line[1])
m_pi2 = m_pi**2

data_corr = np.loadtxt(corr_filename, skiprows=2)
t = data_corr[:,0]
Gt = data_corr[:,1:]


colours = ['red', 'blue']

for ib in range(n_bare):
    ax.plot(t[:], Gt[:,ib], color=colours[ib])


ax.set_title(f'P33 {int(n_bare)}b{int(n_ch)}c, L = {L:.1f} fm, ' + '$m_{\pi}^{2}$' + f' = {m_pi2:.3f} GeV' + '${}^2$')
ax.set_xlabel('$t$ (fm)')
ax.set_ylabel('$G(t)$')

# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

outfile = f'figs/P33_{int(n_bare)}b{int(n_ch)}c_correlation_L{int(L*10)}_mpisq{int(m_pi2*1000)}_fit{fitNum}.{figType}'
if saveFigs:
    savefig(outfile, bbox_inches='tight')
else:
    plt.show()


# ------------------------------------------------------------------
fig, ax = plt.subplots()
cont_filename = f'data/contaminationMulti_fit{fitNum}.out'
data_cont = np.loadtxt(cont_filename, skiprows=2)
t = data_cont[:,0]
Ct = data_cont[:,1:]


for ib in range(n_bare):
    ax.plot(t[:], Ct[:,ib], color=colours[ib])

ax.set_title(f'P33 {int(n_bare)}b{int(n_ch)}c, L = {L:.1f} fm, ' + '$m_{\pi}^{2}$' + f' = {m_pi2:.3f} GeV' + '${}^2$')
ax.set_xlabel('$t$ (fm)')
ax.set_ylabel('$C_{i}(t)$')


# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

outfile = f'figs/P33_{int(n_bare)}b{int(n_ch)}c_contamination_L{int(L*10)}_mpisq{int(m_pi2*1000)}_fit{fitNum}.{figType}'
if saveFigs:
    savefig(outfile, bbox_inches='tight')
else:
    plt.show()
