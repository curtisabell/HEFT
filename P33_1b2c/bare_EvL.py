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
HEFT.printHEFTInfo()

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()

saveFigs = True
m_pi0 = 0.1385

# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])
print(f'Fit {fitNum}')

H_eigs_file = np.loadtxt(f"data/H_eigenvalues_L_fit{fitNum}.out", skiprows=0)
H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_L_fit{fitNum}.out", skiprows=0)
bare_state = np.loadtxt(f"data/bare_state_L_fit{fitNum}.out", skiprows=0)
bare_index = bare_state[:,4:7]
bare_state = bare_state[:,0:4]
Lambda = np.loadtxt(f"data/finiteParams_L_fit{fitNum}.out", skiprows=0)

L = H_eigs_file[:,0]
nL = len(L)
H_eigs = H_eigs_file[:,1:]
nEigs = len(H_eigs[0,:])
ch_labels = H0_eigs_file[0,1:]
H0_eigs = H0_eigs_file[1:,1:]

plotBasis = False

if plotBasis:
    pypl.plot(L, H0_eigs[:,0], color='blue', label='_nolegend_')
    pypl.plot(L, H0_eigs[:,0], color='blue', linestyle='dashed', label='__nolegend__')
    pypl.plot(L, H0_eigs[:,1], color='blue', linestyle='dashed', label='__nolegend__')
pypl.plot(L, H_eigs[:,0], color='black', label='__nolegend__')

nEigsPlot = 20
for i in range(1,nEigsPlot):
    if plotBasis:
        if (i >= 2):
            pypl.plot(L, H0_eigs[:,i], color='blue'
                      , linestyle='dashed', label='_nolegend_')
    pypl.plot(L, H_eigs[:,i], color='black'
              , label='_nolegend_')
j = 0
old_min1 = 0
old_min2 = 0
old_min3 = 0
dddotted = (0, (3, 1, 1, 1))

# this is a stupid method good luck ever remembering
#    how this actually works
while (j < nL):
    if ((bare_index[j,0] != bare_index[j-1,0]) or (j==(nL-1))):
        pypl.plot(L[old_min1:(j)], bare_state[old_min1:(j),1]
                  , color='red', linewidth=4, linestyle='solid'
                  , label='_nolegend')
        old_min1 = j
    if ((bare_index[j,1] != bare_index[j-1,1]) or (j==(nL-1))):
        pypl.plot(L[old_min2:(j)], bare_state[old_min2:(j),2]
                  , color='blue', linewidth=4, linestyle='--', dashes=(5,5)
                  , label='_nolegend')
        old_min2 = j
    if ((bare_index[j,2] != bare_index[j-1,2]) or (j==(nL-1))):
        pypl.plot(L[old_min3:(j)], bare_state[old_min3:(j),3]
                  , color='green', linewidth=4, linestyle='--', dashes=(20,10)
                  , label='_nolegend')
        old_min3 = j
    j = j + 1

# red_patch = mpatches.Patch(color='red', label='1st')
# green_patch = mpatches.Patch(color='green', label='2nd')
# cyan_patch = mpatches.Patch(color='cyan', label='3rd')
# pypl.legend(handles=[red_patch, green_patch, cyan_patch]
#             , labels=['1st most probable', '2nd most probable'
#                       , '3rd most probable'], prop={'size': 20}
#             , loc='upper right')

bare1, = pypl.plot([-1, -1], [-1, -1], color='red',
                   linewidth=4.0, linestyle='solid')
bare2, = pypl.plot([-1, -1], [-1, -1], color='blue',
                   linewidth=4.0, linestyle='--', dashes=(5,5))
bare3, = pypl.plot([-1, -1], [-1, -1], color='green',
                   linewidth=4.0, linestyle='--', dashes=(8,8))

pypl.legend([bare1, bare2, bare3]
            , ['1st most probable', '2nd most probable', '3rd most probable']
            , loc='upper right', fontsize=18)

# pypl.legend(handles=[red_patch]
#             , labels=['Largest $\Delta^{(0)}$ Componant'], prop={'size': 20}
#             , loc='upper right')

# pypl.axis([min(L), max(L), 1.05, 1.5])
pypl.axis([min(L), max(L), 1.05, 1.8])
pypl.ylabel('E (GeV)')
pypl.xlabel('L (fm)')
# pypl.title('$\Lambda = %.2f$ GeV' % Lambda)
# pypl.title('$\Lambda$ = 0.8 GeV')
# pypl.title('Finite-volume energy spectrum')
# pypl.title('1b1c, $v_{\pi N}^{\pi N}$ = 0')
# pypl.title('1b2c')

outFile = 'figs/1b2c_EvL_Lam800MeV_bare.pdf'

if saveFigs:
    print('Writing to file: '+outFile)
    savefig(outFile, bbox_inches='tight')
else:
    pypl.show()
