#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

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

# -----------------------------Read data----------------------------
with open(f'data/latest_finiteProgram_fit{fitNum}.out', 'r') as f:
    finProg = f.readline().strip()
if (finProg=='lqcd'):
    plotVarL = True
elif (finProg=='mpi'):
    plotVarL = False
else:
    print('Error in latest finiteProgram file')
    exit()

if plotVarL:
    H_evec_file = np.loadtxt(f"data/H_eigenvectors_lqcd_fit{fitNum}.out", skiprows=0)
    H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_lqcd_fit{fitNum}.out", skiprows=0)
    H_eigs_file = np.loadtxt(f"data/H_eigenvalues_lqcd_fit{fitNum}.out", skiprows=0)
else:
    H_evec_file = np.loadtxt(f"data/H_eigenvectors_m_pi_fit{fitNum}.out", skiprows=0)
    H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)
    H_eigs_file = np.loadtxt(f"data/H_eigenvalues_m_pi_fit{fitNum}.out", skiprows=0)


m_pis = H0_eigs_file[1:,0]
H_eigs = H_eigs_file[:,1:]
nMpi = int(len(m_pis[:]))
nBasis = int(len(H_evec_file[:,0]) / nMpi)
nEigs = len(H_evec_file[0,1:])
eVects = np.zeros((nMpi, nBasis, nEigs))
ch_numbers = H0_eigs_file[0,1:].astype('int')
n_k = np.count_nonzero(ch_numbers==1)

for i in range(nMpi):
    for j in range(nBasis):
        row = i*nBasis + j
        eVects[i,j,:] = H_evec_file[row,1:]

# nPoints    5
# L          3.0
# sites      32
# a        m_pi     m      dm
# 0.0933   0.169    1.49   0.12
# 0.0950   0.280    1.56   0.05
# 0.0961   0.391    1.71   0.03
# 0.10086  0.515    1.78   0.05
# 0.1022   0.623    1.90   0.04
# a        m_pi     m      dm
# 0.0933   0.169    1.55   0.08
# 0.0950   0.280    1.75   0.09
# 0.0961   0.391    1.77   0.03
# 0.10086  0.515    1.82   0.02
# 0.1022   0.623    1.95   0.02
chosen_a = 0.0961
GeV_fm = 0.1973
t_max = 10.0
nt = 100
t = np.linspace(0.0, t_max, nt) # fm
dt = t[1] - t[0]

chosen_mpi2 = 0.139**2
# chosen_mpi2 = 0.623**2
chosen_index = np.argmin(np.abs(m_pis[:] - chosen_mpi2))

Gt = np.zeros((n_bare, nt))
for i_t in range(nt):
    Gt[0,i_t] = np.sum((eVects[chosen_index,0,:nEigs])**2 \
                         * np.exp(-t[i_t]*H_eigs[chosen_index,:nEigs]/GeV_fm))

# Gt = 0.0
# for i_b in range(n_bare):
#     for i_E in range(nEigs):
#         Gt[i_b,:] = Gt[i_b,:] + (eVects[chosen_index,i_b,i_E])**2 \
#             * np.exp(-H_eigs[chosen_index,i_E]*t[:] / GeV_fm)


# # -----------------------Correlation Functions----------------------
# ax.plot(t, Gt[0,:], color='red')
# ax.set_ylim([0.0, 1.0])
# ax.set_xlim([0.0, t_max])
# ax.set_ylabel('G(t)')
# ax.set_xlabel('$t$ (fm)')


# -------------------------Effective masses-------------------------
# E_eff = 1.0/dt * np.log(Gt[:,:-1] / Gt[:,1:])
# ax.plot(t[1:], E_eff[0,:], color='red')

E_eff = ( Gt[0,:] - eVects[chosen_index,0,0]**2 \
          * exp(-t[:]*H_eigs[chosen_index,0]/GeV_fm)) / Gt[0,:]
ax.plot(t[:], E_eff[:], color='red')

print(E_eff)
ax.set_ylim([1.2, 2.8])
# ax.set_ylim([1.65, 1.8])
ax.set_xlim([0.0, t_max])
ax.set_ylabel('$E_{\mathrm{eff}}$')
ax.set_xlabel('$t$ (fm)')


# ax.set_title('$m_{\pi}^2 = $' + f' {chosen_mpi2:.3f} GeV' + '${}^2$')
ax.set_title('$m_{\pi} = $' + f' {np.sqrt(chosen_mpi2):.3f} GeV')

# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

if saveFigs:
    outfile = f'figs/corr.{figType}'
    savefig(outfile, bbox_inches='tight')
else:
    plt.show()
