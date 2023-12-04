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

# --------------------Read pseudo scattering data-------------------
pseudo = np.loadtxt("dataInf_pseudo.in", skiprows=1)
# If in MeV, convert to GeV
if pseudo[0,0]>=1000:
    pseudo[:,0] = pseudo[:,0] / 1000
Epseudo = pseudo[:,0]
phasePseudo = pseudo[:,1]
phaseErrorPseudo = pseudo[:,2]
SrPseudo = pseudo[:,3]
SrErrorPseudo = pseudo[:,4]
etaPseudo = np.sqrt(ones(len(SrPseudo))-SrPseudo)
etaErrorPseudo = SrErrorPseudo / (2.0*etaPseudo)

# uncert = 0.02
# uncert = 0.05

# phasePseudoUp = phasePseudo + phasePseudo*uncert
# phasePseudoDown = phasePseudo - phasePseudo*uncert

phasePseudoUp = phasePseudo + phaseErrorPseudo
phasePseudoDown = phasePseudo - phaseErrorPseudo

# etaPseudoUp = etaPseudo + etaPseudo*uncert
# etaPseudoDown = etaPseudo - etaPseudo*uncert

etaPseudoUp = etaPseudo + etaErrorPseudo
etaPseudoDown = etaPseudo - etaErrorPseudo


# # -------------------------Plot Phase Shift-------------------------
ax.fill_between(Epseudo, y1=phasePseudoUp, y2=phasePseudoDown, color='cyan')
pypl.plot(Epseudo, phasePseudo, color='blue', linestyle='-', linewidth='1.7', zorder=4.0)
pypl.errorbar(Edata, phaseData, yerr=phaseErrorData, fmt='.',
              ecolor='black', color='black', capsize=2, zorder=5.0)

# eta-N threshold
thres = 0.5479+0.9385 # 1.4864
pypl.plot([thres, thres], [-360, 360]
          , color='black', label='_nolegend_', linestyle='dashed')

# K-Lambda threshold
thres = 0.4937+1.116 # 1.6097
pypl.plot([thres, thres], [-360, 360]
          , color='black', label='_nolegend_', linestyle='dashed')


pypl.axis( [min(Epseudo), max(Epseudo), 0, max(phaseData)*1.2] )
# pypl.axis( [min(E), max(E), 0, 200] )
pypl.xlabel( 'E (GeV)' )
pypl.ylabel( '$\delta_{\pi N}$' )
Lam = 0.8
pypl.title(f'S11 {n_bare}b{n_ch}c pseudo')

# xtickloc = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
# pypl.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
pypl.tick_params(axis='y', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='x', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
pypl.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

pypl.show()


# -------------------------Plot Inelasticity------------------------
pypl.cla()

# Plot thresholds
thres = 0.5479+0.9385 # 1.4864
pypl.plot([thres, thres], [min(etaData)*0.8, 1.05]
          , color='black', label='_nolegend_', linestyle='dashed')

thres = 0.4937+1.116 # 1.6097
pypl.plot([thres, thres], [min(etaData)*0.8, 1.05]
          , color='black', label='_nolegend_', linestyle='dashed')

# Plot data
ax.fill_between(Epseudo, y1=etaPseudoUp, y2=etaPseudoDown, color='cyan')
pypl.plot(Epseudo, etaPseudo, color='blue', linestyle='-', linewidth='1.7', zorder=4.0)
pypl.errorbar(Edata, etaData, yerr=etaErrorData, fmt='.',
              ecolor='black', color='black', capsize=2, zorder=5.0)

pypl.ylabel( '$\eta$' )
pypl.xlabel( '$E$ (GeV)' )
pypl.title(f'S11 {n_bare}b{n_ch}c fit{fitNum}')
pypl.axis( [min(Epseudo), max(Epseudo), min(etaData)*0.8, 1.05] )

# xtickloc = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
# pypl.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
pypl.tick_params(axis='y', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='x', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
pypl.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

pypl.show()
