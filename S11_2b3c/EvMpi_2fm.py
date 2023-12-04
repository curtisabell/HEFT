#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as pypl
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

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()

m_pi0 = 0.1385
# ---------------------------Get fit info---------------------------
# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

# ------------------------Check to save figs------------------------
saveFigs = False
figType = 'pdf'

cmdArgs = sys.argv
if (len(cmdArgs)>1):
    saveArg = str(cmdArgs[1])
    if (saveArg.lower()=='save'):
        saveFigs = True
if (len(cmdArgs)>2):
    typeArg = str(cmdArgs[2])
    if (typeArg.lower()=='pdf'):
        figType = 'pdf'

# --------------------Read mpiFiniteVol.f90 data--------------------
with open(f'data/latest_finiteProgram_fit{fitNum}.out', 'r') as f:
    finProg = f.readline().strip()
if (finProg=='lqcd'):
    plotVarL = True
elif (finProg=='mpi'):
    plotVarL = False
else:
    print('Error in latest finiteProgram file')
    exit()

if plotVarL:
    H_eigs_file = np.loadtxt(f"data/H_eigenvalues_lqcd_fit{fitNum}.out", skiprows=0)
    H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_lqcd_fit{fitNum}.out", skiprows=0)
    bare_data = np.loadtxt(f"data/bare_state_lqcd_fit{fitNum}.out", skiprows=0)
    info = np.loadtxt(f"data/finiteParams_lqcd_fit{fitNum}.out", skiprows=1)
else:
    H_eigs_file = np.loadtxt(f"data/H_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)
    H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)
    bare_data = np.loadtxt(f"data/bare_state_m_pi_fit{fitNum}.out", skiprows=0)
    info = np.loadtxt(f"data/finiteParams_m_pi_fit{fitNum}.out", skiprows=1)

bare_state_all = bare_data[:,0:4]
bare_index_all = bare_data[:,4:7]
n_ch = int(info[0])
n_bare = int(info[1])
L = info[2]
Lam_max = info[3]


doPlotLQCD = True
doPlotXerr = True
lqcd_capsize=4.0
lqcd_linewidth=2.0
langfmt = 'v'
JLabfmt = '^'
# -----------------------Plot the lqcd points-----------------------
if doPlotLQCD:
    with open('slopeFitMasses_2fm.data', 'r') as f:
        nMass = int(f.readline().split()[1])
        f.readline()
        nSites = int(f.readline().split()[1])
        massData = np.zeros((nMass, 4))

        for i in range(2):
            f.readline()
            for nM in range(nMass):
                thisLine = np.array(f.readline().split()).astype(float)
                massData[nM,:] = thisLine

            if i==0:
                # Plot the Lang & Verduci points
                ax.errorbar(massData[0,1]**2, massData[0,2]
                            # , xerr=m_pi2_err if doPlotXerr else 0.0*m_pi2_err
                            , yerr=massData[0,3], fmt=langfmt, ecolor='black'
                            , color='black', capsize=lqcd_capsize
                            , label='_nolegend', zorder=8.0, elinewidth=lqcd_linewidth)

            # Plot the JLab points
            ax.errorbar(massData[1:,1]**2, massData[1:,2]
                        # , xerr=m_pi2_err if doPlotXerr else 0.0*m_pi2_err
                        , yerr=massData[1:,3], fmt=JLabfmt, ecolor='black'
                        , color='black', capsize=lqcd_capsize
                        , label='_nolegend', zorder=8.0, elinewidth=lqcd_linewidth)

    # Plot the third Lang & Verduci point
    # ax.errorbar(np.array([0.266])**2, np.array([1.792])
                # , yerr=[0.029], fmt=langfmt, ecolor='black'
                # , color='black', capsize=lqcd_capsize
                # , label='_nolegend', zorder=8.0, elinewidth=lqcd_linewidth)
# ------------------------------------------------------------------


# -------------------------Plot mpiFin data-------------------------
m_pi2 = H_eigs_file[:,0]
n_m_pi = len(m_pi2)
H_eigs = H_eigs_file[:,1:]
nEigs = len(H_eigs[0,:])
H0_eigs = H0_eigs_file[1:,1:]
ch_num = H0_eigs_file[0,1:].astype(int)

# basisColors = ['gray','blue','red','green']
basisColors = ['gray', 'orange', 'green','darkviolet']
style1 = 'dashed'
style2 = (0, (14, 8))
style3 = (0, (8, 6, 3, 6))
basisStyles = ['solid', style2, style1, style3]

doPlotBasis = False
doPlotEigenvalues = True
for i in range(0,nEigs):
    # Plot free-particle states
    if doPlotBasis:
        if i>1:
            pypl.plot(m_pi2, H0_eigs[:,i], color=basisColors[ch_num[i]]
                      , linestyle=basisStyles[ch_num[i]], label='_nolegend_')

    # Plot eigenvalues
    if doPlotEigenvalues:
        pypl.plot(m_pi2, H_eigs[:,i], color='black'
                  , label='_nolegend_', zorder=2.0)

# Plot physical pion mass
pypl.plot([m_pi0**2, m_pi0**2], [-4.0, 4.0], linestyle='dashed'
          , color='black', label='_nolegend')

bs_colours = ['red', 'blue']

plotBareBasis = False
# Plot bare state first
if (doPlotBasis and plotBareBasis):
    pypl.plot(m_pi2, H0_eigs[:,0], color=bs_colours[0]
              , linestyle='solid', label='_nolegend_', zorder=5.0, linewidth=1.0)
    pypl.plot(m_pi2, H0_eigs[:,1], color=bs_colours[1]
              , linestyle='solid', label='_nolegend_', zorder=5.0, linewidth=1.0)


# ---------------------Bare state line plotting---------------------
doPlotBare = True
if (len(cmdArgs)>3):
    bareArg = str(cmdArgs[3])
    if (bareArg.lower()=='raw'):
        doPlotBare = False


# Separate bare_index_all
nr, nc = bare_state_all.shape
bare_state = np.zeros((n_bare, int(nr/2), int(nc)))
nr, nc = bare_index_all.shape
bare_index = np.zeros((n_bare, int(nr/2), int(nc)))

# rearrange the bare-dominated state data to be easier to plot
bare_state[0,:,:] = bare_state_all[0::2]
bare_state[1,:,:] = bare_state_all[1::2]
bare_index[0,:,:] = bare_index_all[0::2]
bare_index[1,:,:] = bare_index_all[1::2]

longdash = (0, (12, 8))

j = 0
old_min1 = 0
old_min2 = 0
old_min3 = 0
if doPlotBare:
    while (j < n_m_pi):
        if ((bare_index[0,j,0] != bare_index[0,j-1,0]) or (j==(n_m_pi-1))):
            pypl.plot(m_pi2[old_min1:(j)], bare_state[0,old_min1:(j),1]
                      , color=bs_colours[0], linewidth=3, alpha=1.0
                      , zorder=5.0, label='_nolegend')
            old_min1 = j
        if ((bare_index[0,j,1] != bare_index[0,j-1,1]) or (j==(n_m_pi-1))):
            pypl.plot(m_pi2[old_min2:(j)], bare_state[0,old_min2:(j),2]
                      , color=bs_colours[0], linewidth=4, linestyle=longdash
                      , alpha=1.0, zorder=3.0, label='_nolegend')
            old_min2 = j
        j = j + 1

j = 0
old_min1 = 0
old_min2 = 0
old_min3 = 0
if doPlotBare:
    while (j < n_m_pi):
        if ((bare_index[1,j,0] != bare_index[1,j-1,0]) or (j==(n_m_pi-1))):
            pypl.plot(m_pi2[old_min1:(j)], bare_state[1,old_min1:(j),1]
                      , color=bs_colours[1], linewidth=3
                      , alpha=1.0, zorder=4.0, label='_nolegend')
            old_min1 = j
        if ((bare_index[1,j,1] != bare_index[1,j-1,1]) or (j==(n_m_pi-1))):
            pypl.plot(m_pi2[old_min2:(j)], bare_state[1,old_min2:(j),2]
                      , color=bs_colours[1], linewidth=4, linestyle=longdash
                      , alpha=1.0, zorder=4.0, label='_nolegend')
            old_min2 = j
        j = j + 1


# Dummy lines for the legend
if doPlotBare:
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=bs_colours[0]
              , linewidth=3.0, label='Largest $N_{1}$')
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=bs_colours[0]
              , linewidth=4.0, linestyle=longdash
              , label='2nd Largest $N_{1}$')
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=bs_colours[1]
              , linewidth=3.0, label='Largest $N_{2}$')
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=bs_colours[1]
              , linewidth=4.0, linestyle=longdash
              , label='2nd Largest $N_{2}$')

if doPlotBasis:
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=basisColors[1]
              , linewidth=1.0, linestyle=basisStyles[1]
              , label='$\pi N$')
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=basisColors[2]
              , linewidth=1.0, linestyle=basisStyles[2]
              , label='$\eta N$')
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=basisColors[3]
              , linewidth=1.0, linestyle=basisStyles[3]
              , label='$K \Lambda$')

pypl.errorbar([-5.0,-5.0],[-5.0,-5.0], yerr=[0.1,0.1], color='black'
              , ecolor='black', capsize=lqcd_capsize
              , elinewidth=lqcd_linewidth, fmt=langfmt
              , label='Lang \& Verduci')
pypl.errorbar([-5.0,-5.0],[-5.0,-5.0], yerr=[0.1,0.1], color='black'
              , ecolor='black', capsize=lqcd_capsize
              , elinewidth=lqcd_linewidth, fmt=JLabfmt
              , label='HSC')
pypl.legend(loc='lower right', numpoints=1
            , prop={'size': 16}, frameon=False)


# ------------------------------------------------------------------
doTitle = False
if doTitle:
    if plotVarL:
        # pypl.title(f'S11 N{n_bare}b{n_ch}c, L ' + '$\sim$' + f' 2.0 fm ,  Fit {fitNum}')
        pypl.title('L ' + '$\sim$' + f' 2.0 fm')
        # pypl.title('2 bare states')
    else:
        pypl.title(f'S11 N{n_bare}b{n_ch}c, L = {L:.2f} fm,  Fit {fitNum}')

max_E = 2.0
pypl.axis([0.0, 0.4, 1.0, max_E])

pypl.ylabel('E (GeV)')
pypl.xlabel('$m_{\pi}^2 (\\textrm{GeV}^2)$')

# xtickloc = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
# pypl.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
pypl.tick_params(axis='y', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='x', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
pypl.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')


outFileFitNum = False
if saveFigs:
    outFile = f'figs/S11_{n_bare}b{n_ch}c_' \
        + ('EvMpi_' if doPlotEigenvalues else '') \
        + ('bare_' if doPlotBare else '') \
        + ('withBasis_' if doPlotBasis else '') \
        + f'L{int(L)}fm' \
        + (f'_fit{fitNum}' if outFileFitNum else '') \
        + f'.{figType}'
    savefig(outFile, bbox_inches='tight')

else:
    pypl.show()
