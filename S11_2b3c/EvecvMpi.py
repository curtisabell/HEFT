#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys
import os

sys.path.append('../heftCode')
import readHEFTConfig
HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()
n_ch = HEFT.n_ch
n_bare = HEFT.n_bare


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

# -------------------------Get lattice size-------------------------
# with open('HEFTFinite.config', 'r') as f:
#     for line in f.readlines():
#         if (line[0:6]=='L_m_pi'):
#             L_m_pi = float(line.strip().split()[1])
#             Lfm = int(round(L_m_pi))
#             print(f'L = {L_m_pi:.2f} fm')


# -----------------------------Read data----------------------------
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
    H_evec_file = np.loadtxt(f"data/H_eigenvectors_lqcd_fit{fitNum}.out", skiprows=0)
    H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_lqcd_fit{fitNum}.out", skiprows=0)
    info = np.loadtxt(f"data/finiteParams_lqcd_fit{fitNum}.out", skiprows=1)
else:
    H_evec_file = np.loadtxt(f"data/H_eigenvectors_m_pi_fit{fitNum}.out", skiprows=0)
    H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)
    info = np.loadtxt(f"data/finiteParams_m_pi_fit{fitNum}.out", skiprows=1)


m_pis = H0_eigs_file[1:,0]
nMpi = int(len(m_pis[:]))
nBasis = int(len(H_evec_file[:,0]) / nMpi)
nEigs = len(H_evec_file[0,1:])
eVects = np.zeros((nMpi, nBasis, nEigs))
ch_numbers = H0_eigs_file[0,1:].astype('int')
n_k = np.count_nonzero(ch_numbers==1)
L = info[2]
Lfm = int(round(L))


for i in range(nMpi):
    for j in range(nBasis):
        row = i*nBasis + j
        eVects[i,j,:] = H_evec_file[row,1:]


# nEigs = 5 # how many eigenvalues to plot
# H_eigs = H_eigs_file[:,1:nEigs+1]
# H0_eigs = H0_eigs_file[1:,1:nEigs+1]
# ch_numbers = H0_eigs_file[0,1:].astype('int')

# # plot eigenvaluse
# for i in range(0,nEigs):
#     ax.plot(indepParam[:], H0_eigs[:,i], color='blue', linestyle='dashed')
#     ax.plot(indepParam[:], H_eigs[:,i], color='black', linestyle='solid')
# ax.set_xlim(min(indepParam), max(indepParam))
# ax.set_xlabel(indepParamName)

doPlotLQCD = True
lqcd_capsize=4.0
# ------------------Plot first odd parity nucleon-------------------
if doPlotLQCD:
    if Lfm==3:
        massData  = np.loadtxt('mass_N1535.in', skiprows=3)
        m_pi_CSSM = massData[:,0]
        CSSM_fmt = 'o'
    elif Lfm==2:
        lqcd_capsize=4.0
        langfmt = 'v'
        JLabfmt = '^'

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




if (Lfm==2):
    nEigs_plot = 4
elif (Lfm==3):
    nEigs_plot = 6
else:
    nEigs_plot = 6

nBasis_plot = 5
nScattering = nBasis_plot - n_bare

colours = ['red', 'blue', 'orange', 'green','darkviolet']
evecLabels = ['$N_1$', '$N_2$'
              , '$\pi N$', '$\eta N$'
              , '$K\Lambda$']

# lstyles = ['dashed' for i in range(len(nBasis_plot))]
# lstyles[0] = 'solid'
# lstyles[1] = 'solid'
dddotted = (0, (3, 1, 1, 1))
longdash = (0, (12, 8))
longdashdot = (0, (8, 6, 3, 6))
style1 = 'dashed'
style2 = (0, (14, 8))
style3 = (0, (8, 6, 3, 6))
lstyles = ['solid', 'solid', style2, style1, style3]
# lstyles = ['solid', 'solid', 'dashed', longdashdot]
zorders = [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]

for eig in range(nEigs_plot):
    plt.style.use('classic')
    plt.rcParams['axes.formatter.useoffset'] = False
    plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rc('text', usetex=True)
    fig, ax = plt.subplots()

    for basis in range(nBasis_plot):
        if ch_numbers[basis]==0:
            lwidth = 3.0
            y = eVects[:,basis,eig]
        else:
            lwidth = 2.0
            # plot the specific scattering basis state
            # y = eVects[:,basis,eig]
            # plot the sum of all contributions from a channel
            y = np.sum(eVects[:,ch_numbers[:nBasis]==ch_numbers[basis],eig], axis=1)

        ax.plot(m_pis[:], y
                , color=colours[basis], linestyle=lstyles[basis]
                , label=evecLabels[basis], linewidth=lwidth
                , zorder=zorders[basis])

    # plot markers at lQCD masses
    if doPlotLQCD:
        lqcd_height = 0.5
        # s should be set to the errorbar capsize**2
        if Lfm==3:
            ax.scatter(m_pi_CSSM**2, np.ones(shape(m_pi_CSSM))*lqcd_height,
                       marker=CSSM_fmt, s=lqcd_capsize**2, zorder=8.0)
        elif Lfm==2:
            ax.scatter(massData[0,1]**2, np.ones(shape(massData[0,1]))*lqcd_height,
                       marker=langfmt, s=lqcd_capsize**2, zorder=8.0)
            ax.scatter(massData[1:,1]**2, np.ones(shape(massData[1:,1]))*lqcd_height,
                       marker=JLabfmt, s=lqcd_capsize**2, zorder=8.0)


    indepParamName = '$m_{\pi}^2$ (GeV${}^2$)'
    ax.set_xlabel(indepParamName, fontsize=20)
    yaxisLabel = '$| \langle B_i | E_{%i} \\rangle |^2$' % (eig+1)
    # yaxisLabel = '$| \langle E_{%i} | B_i \\rangle |^2$' % (eig+1)
    ax.set_ylabel(yaxisLabel, fontsize=20)
    ax.set_xlim([min(m_pis), max(m_pis)])
    ax.set_ylim([0.0, 1.0])

    # if (eig==5):
    #     legendloc = 'center left'
    # else:
    #     legendloc = 'best'
    legendloc = 'best'
    if (eig==0):
        plt.legend(loc=legendloc, numpoints=1
                   , prop={'size': 20}, frameon=False)

    # xtickloc = []
    # ax.xticks(xtickloc)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
    ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
    ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
    ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

    if saveFigs:
        outFile = f'figs/{n_bare}b{n_ch}c_Evec_{finProg}_L{Lfm}fm_E{eig+1}.{figType}'
        savefig(outFile, bbox_inches='tight')
    else:
        plt.show()
    ax.clear()
