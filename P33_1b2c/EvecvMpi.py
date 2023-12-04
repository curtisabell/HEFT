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

# -------------------------Get lattice size-------------------------
with open('HEFTFinite.config', 'r') as f:
    for line in f.readlines():
        if (line[0:6]=='L_m_pi'):
            L_m_pi = float(line.strip().split()[1])
            print(f'L = {L_m_pi:.2f} fm')


# -----------------------------Read data----------------------------
H_evec_file = np.loadtxt(f"data/H_eigenvectors_m_pi_fit{fitNum}.out", skiprows=0)
H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)

m_pis = H0_eigs_file[1:,0]
nMpi = int(len(m_pis[:]))
nBasis = int(len(H_evec_file[:,0]) / nMpi)
nEigs = len(H_evec_file[0,1:])
eVects = np.zeros((nMpi, nBasis, nEigs))
ch_numbers = H0_eigs_file[0,1:].astype('int')
n_k = np.count_nonzero(ch_numbers==1)

for i in range(nMpi):
    for j in range(nBasis):
        row = i*nBasis + j
        eVects[i,j,:] = H_evec_file[row,1:]


nEigs_plot = 2
nBasis_plot = 3
nScattering = nBasis_plot - n_bare

colours = ['red', 'orange', 'darkviolet']
evecLabels = ['$\Delta_0$', '$\pi N$', '$\pi\Delta$']

dddotted = (0, (3, 1, 1, 1))
longdash = (0, (12, 8))
longdashdot = (0, (8, 6, 3, 6))
style1 = 'dashed'
style2 = (0, (14, 8))
style3 = (0, (8, 6, 3, 6))
style4 = (0, (14, 8, 3, 6))
lstyles = ['solid', style2, style1, style3, style4]
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

    indepParamName = '$m_{\pi}^2$ (GeV${}^2$)'
    ax.set_xlabel(indepParamName)
    yaxisLabel = '$| \langle B_j | E_{%i} \\rangle |^2$' % (eig+1)
    ax.set_ylabel(yaxisLabel)
    ax.set_xlim([min(m_pis), max(m_pis)])
    ax.set_ylim([0.0, 1.0])

    legendloc = 'best'
    if (eig==0):
        plt.legend(loc=legendloc, numpoints=1
                   , prop={'size': 16}, frameon=False)

    # xtickloc = []
    # ax.xticks(xtickloc)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
    ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
    ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
    ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

    if saveFigs:
        outFile = f'figs/{n_bare}b{n_ch}c_Evec_L{int(L_m_pi)}fm_E{eig+1}.{figType}'
        savefig(outFile, bbox_inches='tight')
    else:
        plt.show()
    ax.clear()
