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
indepParam = params[:,indepParamIndex]**2
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

nEigs = 5 # how many eigenvalues to plot
H_eigs = H_eigs_file[:,1:nEigs+1]
H0_eigs = H0_eigs_file[1:,1:nEigs+1]
ch_numbers = H0_eigs_file[0,1:].astype('int')

# plot eigenvalues
for i in range(0,nEigs):
    ax.plot(indepParam[:], H0_eigs[:,i], color='blue', linestyle='dashed')
    ax.plot(indepParam[:], H_eigs[:,i], color='black', linestyle='solid')
ax.set_xlim(min(indepParam), max(indepParam))
# ax.set_xlabel(indepParamName)

ax.set_xlabel('$\Lambda^2$ (GeV${}^2$)')
ax.set_ylabel('E (GeV)')
ax.set_ylim(1.1, 1.6)

ax.set_title('L = 5.00 fm')

# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

savePlots = False
figType = 'pdf'

if savePlots:
    outFile = f'figs/{n_bare}b{n_ch}c_EvParam.{figType}'
    savefig(outFile, bbox_inches='tight')
else:
    plt.show()
