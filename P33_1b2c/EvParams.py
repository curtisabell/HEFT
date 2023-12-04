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

plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, ax = plt.subplots()

# Check if a valid argument
if len(sys.argv)!=3:
    print('Please enter the name of a .params file nad an index for the independent variable as arguments')
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
indepParamName = paramNames[indepParamIndex].replace('_', '') # prevent tex interpretation

# Check that multifitFin.x has been run and outputted the file(s)
H_eigs_fileName = 'data/H_eigenvalues_multifit_' + paramsFileName + '.out'
H0_eigs_fileName = 'data/H0_eigenvalues_multifit_' + paramsFileName + '.out'
if os.path.exists(H_eigs_fileName):
    H_eigs_file = np.loadtxt(H_eigs_fileName)
    H0_eigs_file = np.loadtxt(H0_eigs_fileName)
else:
    print(H_eigs_fileName + ' does not exist')
    print('Please run: ./multifitFin.x ' + paramsFileName)
    sys.exit('Stopping...')

nEigs = 10 # how many eigenvalues to plot
H_eigs = H_eigs_file[:,1:nEigs+1]
H0_eigs = H0_eigs_file[1:,1:nEigs+1]
ch_numbers = H0_eigs_file[0,1:].astype('int')
ch_num = ch_numbers

colours = ['blue', 'orange', 'darkviolet']
labels = ['$\Delta_{0}$', '$\pi N$', '$\pi\Delta$']
dddotted = (0, (3, 1, 1, 1))
longdash = (0, (12, 8))
longdashdot = (0, (8, 6, 3, 6))
style1 = 'dashed'
style2 = (0, (14, 8))
style3 = (0, (8, 6, 3, 6))
style4 = (0, (14, 8, 3, 6))
lstyles = ['dashdot', style1, style2]

# plot eigenvalues
for i in range(0,nEigs):
    lstyle = lstyles[ch_num[i]]
    lcolour = colours[ch_num[i]]
    if i==ch_num[i]:
        llabel = labels[ch_num[i]]
    else:
        llabel = '_nolegend'
    ax.plot(indepParam[:], H0_eigs[:,i], color=lcolour
            , linestyle=lstyle, label=llabel)
    ax.plot(indepParam[:], H_eigs[:,i], color='black', linestyle='solid')
ax.set_xlim(min(indepParam), max(indepParam))

# ax.set_xlabel(indepParamName)
ax.set_xlabel('$\Lambda$ (GeV)')
# ax.set_xlabel('$\Lambda^2$ (GeV${}^2$)')
ax.set_ylabel('E (GeV)')
ax.set_ylim(1.0, 1.8)

ax.set_title('L = 2.99 fm')
# ax.set_title('L = 5.00 fm')

useLegend = True
legendloc = 'lower left'
if useLegend:
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

# savePlots = False
savePlots = True
figType = 'pdf'

if savePlots:
    outFile = f'figs/{n_bare}b{n_ch}c_EvParam.{figType}'
    savefig(outFile, bbox_inches='tight')
else:
    plt.show()
