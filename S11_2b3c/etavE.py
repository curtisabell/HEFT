#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as pypl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()

sys.path.append('../heftCode')
import readHEFTConfig
HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()
n_ch = HEFT.n_ch
n_bare = HEFT.n_bare

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






# Plot thresholds
thres = 0.5479+0.9385 # 1.4864
pypl.plot([thres, thres], [min(etaData)*0.8, 1.05]
          , color='black', label='_nolegend_', linestyle='dashed')

thres = 0.4937+1.116 # 1.6097
pypl.plot([thres, thres], [min(etaData)*0.8, 1.05]
          , color='black', label='_nolegend_', linestyle='dashed')

# Plot data
pypl.errorbar(Edata, etaData, yerr=etaErrorData, fmt='.',
              ecolor='black', color='black', capsize=2)
pypl.plot(E, eta, color='blue', linestyle='-', linewidth='1.7')
pypl.ylabel( '$\eta$' )
pypl.xlabel( '$E$ (GeV)' )
pypl.title(f'S11 {n_bare}b{n_ch}c')
# pypl.title(f'S11 {n_bare}b{n_ch}c fit{fitNum}')
pypl.axis( [min(E), max(E), min(etaData)*0.8, 1.05] )

# xtickloc = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
# pypl.xticks(xtickloc)
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
