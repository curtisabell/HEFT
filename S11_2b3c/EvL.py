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

saveFigs = False
m_pi0 = 0.1385

# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])
print(f'Fit {fitNum}')

# #!/usr/bin/python3
# import matplotlib.pyplot as pypl
# import numpy as np
# from pylab import *
# from matplotlib import rc
# import math

# pypl.style.use('classic')
# pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
# pypl.rc('text', usetex=True)

# fig, ax = pypl.subplots()

H_eigs_file = np.loadtxt(f"data/H_eigenvalues_L_fit{fitNum}.out", skiprows=0)
H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_L_fit{fitNum}.out", skiprows=0)
# bare_state = np.loadtxt("bare_state_L.out", skiprows=0)
# bare_state = bare_state[:,0:3]
Lambda = np.loadtxt(f"data/finiteParams_L_fit{fitNum}.out", skiprows=0)

L = H_eigs_file[:,0]
nL = len(L)
H_eigs = H_eigs_file[:,1:]
nEigs = len(H_eigs[0,:])
ch_labels = H0_eigs_file[0,1:]
H0_eigs = H0_eigs_file[1:,1:]

# pypl.plot(L, H0_eigs[:,0], color='red')
pypl.plot(L, H0_eigs[:,0], color='blue', linestyle='dashed')
pypl.plot(L, H_eigs[:,0], color='black')

nEigsPlot = 20
for i in range(1,nEigsPlot):
    pypl.plot(L, H0_eigs[:,i], color='blue'
                  , linestyle='dashed', label='_nolegend_')
    pypl.plot(L, H_eigs[:,i], color='black'
              , label='_nolegend_')

# pypl.plot([min(L), max(L)], [1.232, 1.232], label='_nolegend_')
pypl.axis([min(L), max(L), 1.1, 1.8])
pypl.ylabel('E (GeV)')
pypl.xlabel('L (fm)')

# xtickloc = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35]
# pypl.xticks(xtickloc)

#  Set the frequency of the axis sub ticks.  Calling AutoMinorLocator() will
#    set the sub ticks automatically, while to have n subintervals
#    call AutoMinorLocator(n).
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

#  Adjust the tick parameters.  Further input parameters can be found at
#  http://matplotlib.org/api/axes_api.html?highlight=tick_params#matplotlib.axes.Axes.tick_params
#
pypl.tick_params(axis='y', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='x', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
pypl.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

# pypl.title('$\Lambda_v$ = %i MeV' % (Lambda*1000))
# pypl.title('$\Lambda_v$ = %i MeV, no bare' % (Lambda*1000))
# pypl.title('$\Lambda_v$ = %i MeV, $\Delta$ = 1450, $\Lambda$ = 800' % (Lambda*1000))
# pypl.title('$\Lambda_v$ = 1000 MeV, no bare')

# savefig('background_finite_%i_MeV.pdf' % (Lambda*1000),bbox_inches='tight')
# savefig('background_finite_bare_800_%i_MeV' % (Lambda*1000),bbox_inches='tight')
# savefig('background_finite_1000_MeV',bbox_inches='tight')

# pypl.title('No background, $\Lambda$ = 1.2 GeV')

# savefig('1b2c_EvL_Lam800MeV.pdf.pdf', bbox_inches='tight')
pypl.show()

# pypl.title('$\Lambda = %.2f$ GeV' % Lambda)
# pypl.title('Finite-volume energy spectrum')
