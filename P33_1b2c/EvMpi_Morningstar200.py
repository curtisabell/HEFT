#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as pypl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()


# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

m_pi0 = 0.139
plotBasis = True
# plotBasis = False

saveFigs = True
# saveFigs = False

# outFileType = 'png'
outFileType = 'pdf'

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

# # Plot lattice qcd data
# pypl.errorbar(m_pi**2, m_Delta, xerr=m_pi2_err, yerr=m_Delta_err,
#               fmt='.', ecolor='black', color='black'
#               , capsize=2, label='_nolegend_', zorder=10.0)

m_pi_MS = 0.2
# E_MS = np.array([6.34, 7.05]) * m_pi_MS
# E_MS_err = np.array([0.06, 0.06]) * m_pi_MS
E_MS = np.array([1.268, 1.410, 1.530, 1.544])
E_MS_err = np.array([0.012, 0.012, 0.012, 0.012])

# Plot lattice qcd data
pypl.errorbar(np.array([m_pi_MS,m_pi_MS,m_pi_MS,m_pi_MS]), E_MS, yerr=E_MS_err,
              fmt='.', ecolor='black', color='black'
              , capsize=2, label='_nolegend_', zorder=10.0)

m_pi2 = H_eigs_file[:,0]
n_m_pi = len(m_pi2)
H_eigs = H_eigs_file[:,1:]
nEigs = len(H_eigs[0,:])
H0_eigs = H0_eigs_file[1:,1:]
ch_num = H0_eigs_file[0,1:].astype(int)


xPlot = [0.0, 1.0]


# # ----------------------Calculate chi2 per dof----------------------
# E1 = H_eigs[:,0]
# E2 = H_eigs[:,0]
# diffs = abs(m_pi2 - m_pi_MS**2)
# minLoc = np.argmin(diffs)
# m_pi2_chi2 = m_pi2[minLoc]
# E_chi2 = np.array([0.0, 0.0])
# E_chi2[0] = E1[minLoc]
# E_chi2[1] = E2[minLoc]



# Plot bare state first
firstLineStyle = 'dashdot'
if (n_bare == 0):
    firstLineStyle = 'dashed'
if (plotBasis):
    yPlot = H0_eigs[-1,0] * np.array([1,1])
    pypl.plot(xPlot, yPlot, color='blue'
              , linestyle=firstLineStyle, label='_nolegend_')

# for i in range(0,nEigs):
nEigs_plot = 7
for i in range(0,nEigs_plot):
    # Plot free-particle states
    if (plotBasis):
        if i>0:
            # pypl.plot(m_pi2, H0_eigs[:,i], color='blue'
            #           , linestyle='dashed', label='_nolegend_')
            yPlot = H0_eigs[-1,i] * np.array([1,1])
            pypl.plot(xPlot, yPlot, color='blue'
                      , linestyle='dashed', label='_nolegend_')
    # Plot eigenvalues
    yPlot = H_eigs[-1,i] * np.array([1, 1])
    # pypl.plot(m_pi2, H_eigs[:,i], color='black'
    #           , label='_nolegend_', zorder=5.0)
    pypl.plot(xPlot, yPlot, color='black'
              , label='_nolegend_', zorder=5.0)

# Bare state line plotting
doPlotBare = False
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

    bare1, = pypl.plot([-1, -1], [-1, -1], color='red'
                       , linestyle='solid', linewidth=4)
    bare2, = pypl.plot([-1, -1], [-1, -1], color='blue'
                       , linestyle='dashed', linewidth=4)
    bare3, = pypl.plot([-1, -1], [-1, -1], color='green'
                       , linestyle=dddotted, linewidth=4)

    # pypl.legend([bare1, bare2, bare3]
    #         , ['1st most probable', '2nd most probable', '3rd most probable']
    #         , loc='lower right', fontsize=20)

x_fac = m_pi_MS * 0.01
max_E = 1.6
max_x = 0.1
pypl.axis([m_pi_MS - x_fac, m_pi_MS + x_fac, 1.15, max_E])


pypl.ylabel('E (GeV)')
pypl.xlabel('$m_{\pi} (\\textrm{GeV})$')

pypl.title('D200, $L$ = 4.16 fm')
# pypl.title('$m_{\pi}$ = 0.2 GeV, $L$ = 4.16 fm, $\chi^2$ = ' + f'{chi2:.1f}')

# xtickloc = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
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

# outFile = f'figs/1b2c_Evmpi_L{L:.2f}fm_Lam{Lam_max:.1f}GeV.pdf'
# outFile = f'figs/1b2c_Evmpi_L{L:.2f}fm_fit{fitNum}.pdf'
# outFile = f'figs/1b2c_Evmpi_L{L:.2f}fm_Lam{Lam_max:.1f}GeV.pdf'

# outFile = f'figs/1b2c_Evmpi_Morningstar_fit{fitNum}_noCol.{outFileType}'
outFile = f'figs/1b2c_Evmpi_D200_m_pi{m_pi_MS:.2f}_fit{fitNum}.{outFileType}'

fig.set_size_inches(5, 8)
pypl.xticks([m_pi_MS])

if saveFigs:
    savefig(outFile, bbox_inches='tight')
else:
    pypl.show()
