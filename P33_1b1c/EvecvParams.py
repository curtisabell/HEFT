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


saveFigs = True
# saveFigs = False
figType = 'pdf'


# Check if a valid argument
if len(sys.argv)!=3:
    print('Please enter the name of a .params file and an index for the independent variable as arguments')
    print('e.g. ./EvecvParams.py renormFits 1')
    sys.exit('Stopping...')
paramsFileName = sys.argv[1]
paramsFileName_full = paramsFileName + '.params'

# Check the params file exists
if os.path.exists(paramsFileName_full):
    print('Using params from ' + paramsFileName_full)
    params = np.loadtxt(paramsFileName_full, skiprows=2)
    nParams = len(params[0,:])
    nFits = len(params[:,0])
else:
    print(paramsFileName_full + ' does not exist')
    sys.exit('Stopping...')

# get names of params
with open('allFits.params', 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    paramLine = f.readline().split()
    paramNames = paramLine[1:nParams+1]

indepParamIndex = int(sys.argv[2]) # which column of params in the independent variable
indepParam = params[:,indepParamIndex]

indepParamName = paramNames[indepParamIndex]
indepParamNameRaw = indepParamName.replace('_', '') # prevent tex interpretation
indepParamName = '$' + indepParamName + '$'


# Check that multifitFin.x has been run and outputted the file(s)
H_evec_fileName = 'data/H_eigenvectors_multifit_' + paramsFileName + '.out'
H0_eigs_fileName = 'data/H0_eigenvalues_multifit_' + paramsFileName + '.out'
if os.path.exists(H_evec_fileName):
    H_evec_file = np.loadtxt(H_evec_fileName)
    H0_eigs_file = np.loadtxt(H0_eigs_fileName)
else:
    print(H_evec_fileName + ' does not exist')
    print('Please run: ./multifitFin.x ' + paramsFileName)
    sys.exit('Stopping...')

nFits = int(len(params[:,0]))
nBasis = int(len(H_evec_file[:,0]) / nFits)
nEigs = len(H_evec_file[0,1:])
eVects = np.zeros((nFits, nBasis, nEigs))
ch_numbers = H0_eigs_file[0,1:].astype('int')

for i in range(nFits):
    for j in range(nBasis):
        row = i*nBasis + j
        eVects[i,j,:] = H_evec_file[row,1:]


# nEigs = 5 # how many eigenvalues to plot
# H_eigs = H_eigs_file[:,1:nEigs+1]
# H0_eigs = H0_eigs_file[1:,1:nEigs+1]
# ch_numbers = H0_eigs_file[0,1:].astype('int')

# # plot eigenvalues
# for i in range(0,nEigs):
#     ax.plot(indepParam[:], H0_eigs[:,i], color='blue', linestyle='dashed')
#     ax.plot(indepParam[:], H_eigs[:,i], color='black', linestyle='solid')
# ax.set_xlim(min(indepParam), max(indepParam))
# ax.set_xlabel(indepParamName)

nEigs_plot = range(4)
nBasis_plot = range(5)

# colours = ['blue','red','orange','yellow','green']
# evecLabels = ['$\Delta^{(0)}$', '$\pi N(k=1)$', '$\pi N(k=2)$'
#               , '$\pi N(k=3)$', '$\pi N(k=4)$', '$\pi N(k=5)$']
colours = ['red', 'orange', 'green', 'blue', 'darkviolet']
evecLabels = ['$\Delta_0$', '$\pi N(k_1)$', '$\pi N(k_2)$'
              , '$\pi N(k_3)$', '$\pi N(k_{\geq 4})$' ]

# lstyles = ['dashed' for i in range(len(nBasis_plot))]
# lstyles[0] = 'solid'

# lstyles = ['dashed' for i in range(len(nBasis_plot))]
# lstyles[0] = 'solid'
# lstyles[1] = 'solid'
dddotted = (0, (3, 1, 1, 1))
longdash = (0, (12, 8))
longdashdot = (0, (8, 6, 3, 6))
style1 = 'dashed'
style2 = (0, (14, 8))
style3 = (0, (8, 6, 3, 6))
style4 = (0, (14, 8, 3, 6))
lstyles = ['solid', style2, style1, style3, style4]
zorders = [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]
# lstyles = ['solid', 'solid', 'dashed', longdashdot]

# custom settings
indepParamName = '$\Lambda^2$ (GeV${}^2$)'
indepParam = indepParam**2




for eig in nEigs_plot:
    plt.style.use('classic')
    plt.rcParams['axes.formatter.useoffset'] = False
    plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rc('text', usetex=True)
    fig, ax = plt.subplots()

    for basis in nBasis_plot:
        if ch_numbers[basis]==0:
            lwidth = 3.0
        else:
            lwidth = 2.0

        y = eVects[:,basis,eig]
        # plot the remaining ones for the last one
        if basis==(max(nBasis_plot)):
            y = np.sum(eVects[:,basis:,eig], axis=1)


        # ax.plot(indepParam[:], eVects[:,basis,eig]
        ax.plot(indepParam[:], y
                , color=colours[basis], linestyle=lstyles[basis]
                , label=evecLabels[basis], linewidth=lwidth
                , zorder=zorders[basis])


    ax.set_xlabel(indepParamName)
    yaxisLabel = '$| \langle B_j | E_{%i} \\rangle |^2$' % (eig+1)
    ax.set_ylabel(yaxisLabel)
    ax.set_xlim([min(indepParam), max(indepParam)])
    ax.set_ylim([0.0, 1.0])

    if eig==0:
        plt.legend(loc='best', numpoints=1
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
        outFile = f'figs/{n_bare}b{n_ch}c_EvecvParam_E{eig+1}.{figType}'
        savefig(outFile, bbox_inches='tight')
    else:
        plt.show()
    ax.clear()
