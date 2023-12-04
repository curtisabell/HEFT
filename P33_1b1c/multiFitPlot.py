#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as pypl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys
import subprocess

sys.path.append('../heftCode')
import readHEFTConfig

HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()
n_ch = HEFT.n_ch
n_bare = HEFT.n_bare
HEFT.printHEFTInfo()

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()

saveFigs = True
# saveFigs = False

m_pi0 = 0.1385

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

# Need to get the details of allFits.params
with open('allFits.params', 'r') as f:
    paramLines = f.readlines()
nChoices = len(paramLines) - 5
choice = paramLines[1].split()

# Compile fitBare.x
fitBareSlopes = True
if fitBareSlopes:
    proc = subprocess.call('make fitBare.x', shell=True)
proc = subprocess.call('make mpiFin.x', shell=True)


# fitRange = range(1, nChoices+1)
# fitRange = range(2, 10)
# fitRange = [2, 1, 8, 10, 11] # Dipoles
fitRange = [19, 27, 35, 54, 58] # Gaussians
nFits = len(fitRange)

colours = ['blue', 'green', 'red', 'cyan', 'magenta']

with open('HEFTFinite.config', 'r') as f:
    finiteLines = f.readlines()
pointsLine = finiteLines[13].split()
nPoints = int(pointsLine[1])
E_ground = np.zeros((nFits, nPoints))
E_second = np.zeros((nFits, nPoints))
E_third  = np.zeros((nFits, nPoints))

alphaOrder = 2
alpha = np.zeros(alphaOrder)

allPlots = []
fitLeg = []
for iFit, fitNum in enumerate(fitRange):

    # ---------------------Get the data for this fit--------------------
    newChoice = '    ' + choice[0] + f'       {fitNum}\n'
    paramLines[1] = newChoice

    # Write new param choice to file
    with open('allFits.params', 'w') as f:
        f.writelines(paramLines)

    # Get details of this param set
    thisParamLine = paramLines[fitNum+3].split()
    m_bare = float(thisParamLine[1])
    Lam = float(thisParamLine[3])

    # Get this param set's bare mass slope
    if fitBareSlopes:
        proc = subprocess.call('./fitBare.x', shell=True)

    with open(f'data/bare_mass_{alphaOrder}slope_fit{fitNum}.fit', 'r') as f:
        f.readline()
        # Need to gen for more bare masses
        slopeLine = f.readline().split()
        alpha = np.array(slopeLine[1:alphaOrder+1]).astype(float)

    proc = subprocess.call('./mpiFin.x', shell=True)

    H_eigs_file = np.loadtxt(f"data/H_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)
    if (fitNum == fitRange[0]):
        m_pi2 = H_eigs_file[:,0]
    E_ground[iFit, :] = H_eigs_file[:,1]
    E_second[iFit, :] = H_eigs_file[:, 2]
    E_third[iFit, :]  = H_eigs_file[:, 3]

    # ---------------------------Plot this fit--------------------------
    thisPlot, = pypl.plot(m_pi2[:], E_ground[iFit,:],
                           label='_nolegend_', zorder=5.0)

    plotHigherLevels = True
    if plotHigherLevels:
        pypl.plot(m_pi2[:], E_second[iFit,:], color=colours[iFit],
                  label='_nolegend_', zorder=5.0)
        pypl.plot(m_pi2[:], E_third[iFit,:], color=colours[iFit],
                  label='_nolegend_', zorder=5.0)
    allPlots.append(thisPlot)

    # Detailed legend
    # legString = '$m_{B}=$' + f'{m_bare:.3f}, ' + '$\Lambda=$' \
    #     + f'{Lam:.1f}, '
    # for i,alph in enumerate(alpha):
    #     legString = legString + '$\\alpha$' + f'{(i+1)*2}={alph:.2f}   '

    # Simple legend
    legString = '$\Lambda=$' + f' {Lam:.1f} GeV'
    fitLeg.append(legString)

    # if (fitNum==1):
    #     fitLeg.append('$m_{B}=$' + f'{m_bare:.3f}, ' + '$\\alpha=$' + f'{alpha:.2f}')
    # else:
    #     fitLeg.append('$m_{B}=$' + f'{m_bare:.3f}, ' + '$\Lambda=$'
    #                   + f'{Lam:.1f}, ' + '$\\alpha=$' + f'{alpha:.2f}')



# Plot physical pion mass
pypl.plot([m_pi0**2, m_pi0**2], [-4.0, 4.0], linestyle='dashed',
          color='black', label='_nolegend')

pypl.ylabel('E (GeV)')
pypl.xlabel('$m_{\pi}^2 (\\textrm{GeV}^2)$')

min_E = 1.15
max_E = 1.8
pypl.axis([0.0, max(m_pi2), min_E, max_E])

# pypl.title('Gaussian Regulator')

# pypl.title(f'1b1c')

# Resize plot and chuck legend to the right
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# pypl.legend(allPlots, fitLeg, loc='center left',
#             bbox_to_anchor=(1, 0.5), fontsize = 14)

# Add legend onto plot
pypl.legend(allPlots, fitLeg, numpoints=1, loc='lower right'
            , prop={'size': 16}, frameon=False)


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

# fileType = 'png'
fileType = 'pdf'

# outFile = f'figs/1b1c_Evmpi_multifit_{alphaOrder}slope.{fileType}'
outFile = f'figs/1b1c_Evmpi_multifit_{alphaOrder}slope_Gaussian.{fileType}'

if saveFigs:
    savefig(outFile, bbox_inches='tight')
else:
    pypl.show()
