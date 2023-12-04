#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as pypl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math

savePlots = False

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)

#  We want one figure and some axes.
#
fig, ax = pypl.subplots()

# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])


data = np.loadtxt("dataInf_pseudo.in", skiprows=1)
# If in MeV, convert to GeV
if data[0,0]>=1000:
    data[:,0] = data[:,0] / 1000
Edata = data[:,0]
phaseData = data[:,1]
phaseErrorData = data[:,2]
SrData = data[:,3]
SrErrorData = data[:,4]

uncert = 0.03
phaseErrorData = phaseData * uncert

SrErrorData = SrData * uncert
for dsr in SrErrorData:
    dsr = max(dsr, 0.01)


etaData = np.sqrt(ones(len(SrData))-SrData)
# etaErrorData = SrErrorData / (2.0*etaData)
etaErrorData = etaData * uncert


output = np.loadtxt(f"data/scattering_fit{fitNum}.out", skiprows=1)
E = output[:,0]*1000 / 1000
phase = output[:,1]
eta = output[:,3]

thres = 1.232+0.1385

chi2 = 13.23
# chi2dof = 0.61
chi2dof = chi2 / (56 - 6)

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
Lam = 0.8

pypl.title(f'{int(uncert*100)}\% uncertainty, ' + '$\chi^2_{dof} = $' + f' {chi2dof:.3f}')
# pypl.title('$\Lambda = %.2f$ GeV' % Lam)
# pypl.title('0 bare, 7.85 GeV background')
# savefig('noBare_phase_%iMeV' % int(Lam*1000), bbox_inches='tight')
# savefig('phaseShift_1bare7850MeVbackground' ,bbox_inches='tight')

# savefig('1b2c_phaseShift_Lam800MeV.pdf', bbox_inches='tight')
# savefig('1b2c_phaseShift_LamFree.pdf', bbox_inches='tight')
pypl.show()


# -------------------------Plot Inelasticity------------------------
pypl.cla()
pypl.plot([thres, thres], [0.7, 1.1], color='black', label='_nolegend_', linestyle='dashed')
pypl.errorbar(Edata, etaData, yerr=etaErrorData, fmt='.',
              ecolor='black', color='black', capsize=2)
pypl.plot(E, eta, color='red', linestyle='-', linewidth='1.7')
pypl.ylabel( '$\eta$' )
# pypl.title( '$\pi N$ Inelasticity' )
pypl.xlabel( '$E$ (GeV)' )
pypl.axis( [min(E), max(E), 0.8, 1.05] )
pypl.title(f'{int(uncert*100)}\% uncertainty, ' + '$\chi^2_{dof} = $' + f' {chi2dof:.3f}')
# savefig('1b2c_inelasticity_LamFree.pdf', bbox_inches='tight')
pypl.show()
