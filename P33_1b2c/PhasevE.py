#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as pypl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import os
import sys


# ------------------------------------------------------------------
sys.path.append('../src')
import readHEFTConfig
HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()
n_ch = HEFT.n_ch
n_bare = HEFT.n_bare


# -----------------------Stuff for saving figs----------------------
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


# --------------------------Initialise fig--------------------------
pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()


# ---------------------get fit number from file---------------------
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])


# -----------------------------Read data----------------------------
data = np.loadtxt("dataInf.in", skiprows=1)
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

output = np.loadtxt(f"data/scattering_fit{fitNum}.out", skiprows=1)
E = output[:,0]
phase = output[:,1]
eta = output[:,3]

thres = 1.232+0.1385

# -------------------------Plot Phase Shift-------------------------
pypl.plot(E, phase, color='red', linestyle='-', linewidth='1.7')
pypl.errorbar(Edata, phaseData, yerr=phaseErrorData, fmt='.',
              ecolor='black', color='black', capsize=2, zorder=5.0)
pypl.plot([min(E), max(E)], [90, 90], color='black'
          ,linestyle='dashed', label='_nolegend_')
pypl.plot([thres, thres], [-360, 360], color='black', label='_nolegend_', linestyle='dashed')
pypl.axis( [min(E), max(E), 0, 200] )
pypl.xlabel( 'E (GeV)' )
pypl.ylabel( '$\delta_{\pi N}$' )
Lam = 1.2
pypl.title('$\Lambda = %.1f$ GeV' % Lam)


ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
pypl.tick_params(axis='y', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='x', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
pypl.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

outfile = f'figs/{n_bare}b{n_ch}c_phaseShift_fit{fitNum}.{figType}'
if saveFigs:
    savefig(outfile, bbox_inches='tight')
else:
    pypl.show()





# -------------------------Plot Inelasticity------------------------
ax.clear()
pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()


pypl.plot([thres, thres], [0.7, 1.1], color='black', label='_nolegend_', linestyle='dashed')
pypl.errorbar(Edata, etaData, yerr=etaErrorData, fmt='.',
              ecolor='black', color='black', capsize=2)
pypl.plot(E, eta, color='red', linestyle='-', linewidth='1.7')
pypl.ylabel( '$\eta$' )
# pypl.title( '$\pi N$ Inelasticity' )
pypl.xlabel( '$E$ (GeV)' )
pypl.axis( [min(E), max(E), 0.8, 1.05] )

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
pypl.tick_params(axis='y', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='x', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
pypl.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

outfile = f'figs/{int(n_bare)}b{int(n_ch)}c_eta_fit{fitNum}.{figType}'
if saveFigs:
    savefig(outfile, bbox_inches='tight')
else:
    pypl.show()
