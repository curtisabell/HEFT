#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as pypl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

sys.path.append('../src')
import readHEFTConfig

HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()
HEFT.printHEFTInfo()

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()

# saveFigs = True
saveFigs = False

m_pi0 = 0.1385

# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])
print(f'Fit {fitNum}')

H_eigs_file = np.loadtxt(f"data/H_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)
H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)
bare_data = np.loadtxt(f"data/bare_state_m_pi_fit{fitNum}.out", skiprows=0)
bare_index = bare_data[:, 4:7]
bare_state = bare_data[:, 0:4]
info = np.loadtxt(f"data/finiteParams_m_pi_fit{fitNum}.out", skiprows=1)
n_ch = int(info[0])
n_bare = int(info[1])
L = info[2]
Lam_max = info[3]


m_pi2 = H_eigs_file[:,0]
n_m_pi = len(m_pi2)
H_eigs = H_eigs_file[:,1:]
nEigs = len(H_eigs[0,:])
H0_eigs = H0_eigs_file[1:,1:]
ch_num = H0_eigs_file[0,1:].astype(int)


plotBasis = True
# Plot bare state first
firstLineStyle = 'dashdot'
if (n_bare == 0):
    firstLineStyle = 'dashed'
if (plotBasis):
    pypl.plot(m_pi2, H0_eigs[:,0], color='blue',
              linestyle=firstLineStyle, label='_nolegend_')

nEigs_plot = 6
for i in range(0,nEigs_plot):
    # Plot free-particle states
    if i>0 and plotBasis:
        pypl.plot(m_pi2, H0_eigs[:,i], color='blue'
                  , linestyle='dashed', label='_nolegend_')
    # Plot eigenvalues
    pypl.plot(m_pi2, H_eigs[:,i], color='black'
              , label='_nolegend_', zorder=5.0)

# Plot physical pion mass
pypl.plot([m_pi0**2, m_pi0**2], [-4.0, 4.0], linestyle='dashed',
          color='black')

# Bare state line plotting
doPlotBare = True
plotBareLegend = False
if (n_bare == 0):
    doPlotBare = False

j = 0
old_min1 = 0
old_min2 = 0
old_min3 = 0
dddotted = (0, (3, 1, 1, 1))
if doPlotBare:
    while (j < n_m_pi):
        if ((bare_index[j, 0] != bare_index[j-1, 0]) or (j == (n_m_pi-1))):
            pypl.plot(m_pi2[old_min1:j], bare_state[old_min1:j, 1],
                      color='red', linewidth=4, linestyle='solid')
            old_min1 = j
        if ((bare_index[j,1] != bare_index[j-1, 1]) or (j == (n_m_pi-1))):
            pypl.plot(m_pi2[old_min2:(j)], bare_state[old_min2:(j),2]
                      , color='blue', linewidth=4, linestyle='dashed')
            old_min2 = j
        if ((bare_index[j,2] != bare_index[j-1,2]) or (j==(n_m_pi-1))):
            pypl.plot(m_pi2[old_min3:(j)], bare_state[old_min3:(j),3]
                      , color='green', linewidth=4, linestyle=dddotted)
            old_min3 = j
        j = j + 1

    bare1, = pypl.plot([-1, -1], [-1, -1], color='red', linestyle='solid')
    bare2, = pypl.plot([-1, -1], [-1, -1], color='blue', linestyle='dashed')
    bare3, = pypl.plot([-1, -1], [-1, -1], color='green', linestyle=dddotted)
    if plotBareLegend:
        pypl.legend([bare1, bare2, bare3]
                , ['1st most probable', '2nd most probable', '3rd most probable']
                , loc='lower right', fontsize=20)

# Lam_max = 0.8

# pypl.title('$\Delta$ %db1c, L = %.2f fm, $\Lambda$ = %.2f GeV'
#            % (n_bare, L, Lam_max))

# pypl.title('$\Delta$ %db1c, L = %.2f fm, $\Lambda$ = %.2f GeV'
#            % (n_bare, L, Lam_max))
# pypl.title('$\Delta$ 1b1c, L = %.2f fm' % L)
# pypl.title('$\Delta$ 1b1c, L = %.2f fm, $\Lambda$ = 0.8 GeV' % L)
# pypl.title('2 Channel $\Delta$ spectrum,   L = %.2f fm' % L)

# pypl.title(f'1b2c, $\Lambda$ free')


# pypl.title(f'1b2c, $\Lambda = $ {Lam_max:.1f} GeV, $\\alpha = $ 0.71')


min_E = 1.15
max_E = 1.8
pypl.axis([0.0, max(m_pi2), min_E, max_E])
# pypl.axis([0.0, max(m_pi2[H_eigs[:,0]<=max_E]), 1.15, max_E])
# pypl.axis([0.0, 0.4, 1.15, max_E])

pypl.ylabel('E (GeV)')
pypl.xlabel('$m_{\pi}^2 (\\textrm{GeV}^2)$')

# xtickloc = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35]
# pypl.xticks(xtickloc)

#  Set the frequency of the axis sub ticks.  Calling AutoMinorLocator() will
#    set the sub ticks automatically, while to have n subintervals
#    call AutoMinorLocator(n).
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

#  Adjust the tick parameters.  Further input parameters can be found at
#  http://matplotlib.org/api/axes_api.html?highlight=tick_params#matplotlib.axes.Axes.tick_params
#
pypl.tick_params(axis='y', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='x', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
pypl.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

fileType = 'png'
# fileType = 'pdf'

outFile = f'figs/1b2c_Evmpi_L{L:.2f}fm_Lam{Lam_max:.1f}GeV_fit{fitNum}.{fileType}'

if saveFigs:
    savefig(outFile, bbox_inches='tight')
else:
    pypl.show()
