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


m_pi0 = 0.1385


# ---------------------get fit number from file---------------------
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
    fitNum = int(nFitLine[-3:])
    # also get the Lambda value for this fit
    f.readline()
    f.readline()
    for iLine in range(1,fitNum+1):
        line = f.readline()
        if iLine==fitNum:
            lineArr = line.split()
Lam = float(lineArr[3])
print(f'Fit {fitNum}')

# ------------------------get bare mass slope-----------------------
# get the mass slope order
with open('HEFTFinite.config', 'r') as f:
    for line in f.readlines():
        lineArr = line.strip().split()
        if shape(lineArr)!=(0,):
            if lineArr[0]== 'slopeOrder':
                slopeOrder = int(line.split()[1])

# get the actual slope(s)
massSlopes = np.zeros(slopeOrder)
with open(f'data/bare_mass_'
          + (f'{slopeOrder}' if slopeOrder>1 else f'{slopeOrder}')
          + f'slope_fit{fitNum}.fit', 'r') as f:
    f.readline()
    line = f.readline().strip().split()
    massSlopes = np.array(line[1:]).astype('float')


# ------------------------------------------------------------------
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

massData = np.loadtxt('DeltaMass.in', skiprows=3)
hbarc       = 0.1973 # GeVfm
m_pi        = massData[:,0]
m_pi_err    = massData[:,1]
m_pi2_err   = m_pi_err * 2.0 * m_pi
m_Delta     = massData[:,2]
m_Delta_err = massData[:,3]

# Plot lattice qcd data
pypl.errorbar(m_pi**2, m_Delta, xerr=m_pi2_err, yerr=m_Delta_err,
              fmt='.', ecolor='black', color='black'
              , capsize=2, label='_nolegend_', zorder=10.0)

m_pi2 = H_eigs_file[:,0]
n_m_pi = len(m_pi2)
H_eigs = H_eigs_file[:,1:]
nEigs = len(H_eigs[0,:])
H0_eigs = H0_eigs_file[1:,1:]
ch_num = H0_eigs_file[0,1:].astype(int)



# ----------------------Calculate chi2 per dof----------------------
dof = int(len(m_pi) - 1)
E1 = H_eigs[:,0]
E_chi2 = np.zeros(len(m_pi))
dE_chi2 = np.zeros(len(m_pi))
m_pi2_chi2 = np.zeros(len(m_pi))

for i,m2 in enumerate(m_pi[:]**2):
    diffs = abs(m_pi2-m2)
    minLoc = np.argmin(diffs)
    m_pi2_chi2[i] = m_pi2[minLoc]
    E_chi2[i] = E1[minLoc]
    dE_chi2[i] = (E1[minLoc+1] - E1[minLoc]) / (m_pi2[minLoc+1] - m_pi2[minLoc])

    # print(f'{i+1}   {m_Delta[i]}   {m_Delta_err[i]}    {E_chi2[i]}')


chi2 = sum( (E_chi2[:] - m_Delta[:])**2 / (m_Delta_err[:]**2) )
chi2_dof = chi2 / dof
print()
print(f'chi2 = {chi2}')
print(f'chi2/dof = {chi2_dof}')


plotBasis = True
# Plot bare state first
firstLineStyle = 'dashdot'
if (n_bare == 0):
    firstLineStyle = 'dashed'
if (plotBasis):
    pypl.plot(m_pi2, H0_eigs[:,0], color='blue'
              , linestyle=firstLineStyle, label='_nolegend_')

nEigs_plot = 6
for i in range(0,nEigs_plot):
    # Plot free-particle states
    if (i>0) and (plotBasis):
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
                      , color='blue', linewidth=4, linestyle='--', dashes=(5,5))
            old_min2 = j
        if ((bare_index[j,2] != bare_index[j-1,2]) or (j==(n_m_pi-1))):
            pypl.plot(m_pi2[old_min3:(j)], bare_state[old_min3:(j),3]
                      , color='green', linewidth=4, linestyle='--', dashes=(20,10))
            old_min3 = j
        j = j + 1

    bare1, = pypl.plot([-1, -1], [-1, -1], color='red', linestyle='solid')
    bare2, = pypl.plot([-1, -1], [-1, -1], color='blue', linestyle='--', dashes=(5,5))
    bare3, = pypl.plot([-1, -1], [-1, -1], color='green', linestyle='--', dashes=(8,8))
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

# pypl.title(f'1b1c, $\Lambda$ free')
# pypl.title(f'1b1c, $\Lambda = $ {Lam_max:.1f} GeV, $\\alpha = $ 0.77')
alphaDelta = massSlopes[0]

if slopeOrder==1:
    pypl.text(0.26, 1.24, f'$\Lambda = $'+f' {Lam:.1f} GeV', fontsize=18)
    pypl.text(0.26, 1.19, '$\\alpha_{\Delta_{0}}$' \
              + f' = {alphaDelta:.3f} ' + 'GeV ${}^{-1}$', fontsize=18)
elif slopeOrder==2:
    pypl.text(0.25, 1.29, f'$\Lambda = $'+f' {Lam:.1f} GeV', fontsize=18)
    pypl.text(0.25, 1.24, '$\\alpha_{\Delta_{0},2}$' \
              + f' = {massSlopes[0]:.3f} ' + 'GeV ${}^{-1}$', fontsize=18)
    pypl.text(0.25, 1.19, '$\\alpha_{\Delta_{0},4}$' \
              + f' = {massSlopes[1]:.3f} ' + 'GeV ${}^{-3}$', fontsize=18)
else:
    pypl.text(0.26, 1.24, f'$\Lambda = $'+f' {Lam:.1f} GeV', fontsize=18)

# pypl.text(0.28, 1.18, f'$\chi^2$/d.o.f. = {chi2_dof:.1f}', fontsize=18)

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

fileType = 'pdf'
saveFitNum = False
outFile = f'figs/1b1c_EvMpi_' \
    + f'{slopeOrder}slope_L{L:.2f}fm_' \
    + f'Lam{Lam_max:.1f}GeV' \
    + (f'_fit{fitNum}' if saveFitNum else '') \
    + f'.{fileType}'

if saveFigs:
    savefig(outFile, bbox_inches='tight')
else:
    pypl.show()
