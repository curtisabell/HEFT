#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, ax = plt.subplots()

sys.path.append('../src')
import readHEFTConfig
HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()
n_ch = HEFT.n_ch
n_bare = HEFT.n_bare
projName = HEFT.projectName
print(projName)

# -----------------------Stuff for saving figs----------------------
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

# ---------------------Get fit number from file---------------------
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

# -----------------------Read scattering data-----------------------
plotPseudo = False
if plotPseudo:
    data = np.loadtxt("dataInf_pseudo.in", skiprows=1)
else:
    data = np.loadtxt("dataInf_orig.in", skiprows=1)

# If in MeV, convert to GeV
if data[0,0]>=1000:
    data[:,0] = data[:,0] / 1000
Edata = data[:,0]
phaseData = data[:,1]
phaseErrorData = data[:,2]
SrData = data[:,3]
SrErrorData = data[:,4]
etaData = np.sqrt(ones(len(SrData))-SrData)
etaErrorData = SrErrorData / (2.0*etaData)

# ------------------Read data from infiniteVol.f90------------------
output = np.loadtxt(f"data/scattering_fit{fitNum}.out", skiprows=1)
E = output[:,0]*1000 / 1000
phase = output[:,1]
eta = output[:,n_ch+1]

# -------------------------Plot Phase Shift-------------------------
plt.plot(E, phase, color='blue', linestyle='-', linewidth='1.7')
plt.errorbar(Edata, phaseData, yerr=phaseErrorData, fmt='.',
              ecolor='black', color='black', capsize=2, zorder=5.0)

# # eta-N threshold
# thres = 0.5479+0.9385
# plt.plot([thres, thres], [-360, 360]
#           , color='black', label='_nolegend_', linestyle='dashed')

# # K-Lambda threshold
# thres = 0.4937+1.116
# plt.plot([thres, thres], [-360, 360]
#           , color='black', label='_nolegend_', linestyle='dashed')

plt.axis( [min(E), max(E), min([0, min(phase)]), max(phase)*1.05] )
plt.xlabel( 'E (GeV)' )
plt.ylabel( '$\delta_{\pi N}$' )
Lam = 0.8
plt.title(projName.replace('_',' '))

# xtickloc = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
# plt.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(axis='y', which='major', width=1, length=10, color='black')
plt.tick_params(axis='x', which='major', width=1, length=10, color='black')
plt.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
plt.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

outfile = f'figs/{n_bare}b{n_ch}c_phaseShift_fit{fitNum}.{figType}'
if saveFigs:
    savefig(outfile, bbox_inches='tight')
else:
    # print()
    plt.show()

if (n_ch>1):
    # -------------------------Plot Inelasticity------------------------
    ax.clear()
    plt.style.use('classic')
    plt.rcParams['axes.formatter.useoffset'] = False
    plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rc('text', usetex=True)
    fig, ax = plt.subplots()

    # # Plot thresholds
    # thres = 0.5479+0.9385 # 1.4864
    # plt.plot([thres, thres], [min(etaData)*0.8, 1.05]
    #           , color='black', label='_nolegend_', linestyle='dashed')

    # thres = 0.4937+1.116 # 1.6097
    # plt.plot([thres, thres], [min(etaData)*0.8, 1.05]
    #           , color='black', label='_nolegend_', linestyle='dashed')

    # Plot data
    plt.errorbar(Edata, etaData, yerr=etaErrorData, fmt='.',
                  ecolor='black', color='black', capsize=2)
    plt.plot(E, eta, color='blue', linestyle='-', linewidth='1.7')
    plt.ylabel( '$\eta$' )
    plt.xlabel( '$E$ (GeV)' )
    plt.title(projName.replace('_',' '))
    plt.axis( [min(E), max(E), min(etaData)*0.8, 1.05] )

    # xtickloc = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
    # plt.xticks(xtickloc)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.tick_params(axis='y', which='major', width=1, length=10, color='black')
    plt.tick_params(axis='x', which='major', width=1, length=10, color='black')
    plt.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
    plt.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

    outfile = f'figs/{int(n_bare)}b{int(n_ch)}c_eta_fit{fitNum}.{figType}'
    if saveFigs:
        savefig(outfile, bbox_inches='tight')
    else:
        plt.show()
