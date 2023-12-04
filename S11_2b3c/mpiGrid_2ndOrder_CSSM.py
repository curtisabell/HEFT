#!/usr/bin/python3
import matplotlib.pyplot as pypl
import matplotlib.patches as mpatches
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from pylab import *
from matplotlib import rc
import math
import io
import subprocess
import sys


# In mpiFiniteVol.f90, change readSlopesFromFile to .true.
# Make sure bareMassSlopes.in exists

varL = True
if varL:
    programName = 'lqcdFin.x'
    dataVar = 'lqcd'
else:
    programName = 'mpiFin.x'
    dataVar = 'm_pi'


sys.path.append('../heftCode')
import readHEFTConfig
HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()
n_ch = HEFT.n_ch
n_bare = HEFT.n_bare

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)

m_pi0 = 0.1385
# ---------------------------Get fit info---------------------------
# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

# ------------------Setup mass slopes to loop over------------------
firstOrder = [0.96, 0.75]

slopeLimits = [-0.6, 0.0]
nSlopes = 12

slopes = np.linspace(slopeLimits[0], slopeLimits[1], nSlopes)
slopePairs = []
for xx in range(nSlopes):
    for yy in range(nSlopes):
        slopePairs.append([slopes[xx], slopes[yy]])


massData_N1  = np.loadtxt('mass_N1535.in', skiprows=3)
m_pi_N1      = massData_N1[:,0]
m_pi_N1_err  = massData_N1[:,1]
m_N1         = massData_N1[:,2]
m_N1_err     = massData_N1[:,3]
m_pi2_N1_err = m_pi_N1_err * 2.0 * m_pi_N1

massData_N2  = np.loadtxt('mass_N1650.in', skiprows=3)
m_pi_N2      = massData_N2[:,0]
m_pi_N2_err  = massData_N2[:,1]
m_N2         = massData_N2[:,2]
m_N2_err     = massData_N2[:,3]
m_pi2_N2_err = m_pi_N2_err * 2.0 * m_pi_N2

counter = 0
for bareSlopes in slopePairs:
    counter += 1
    print(f'{counter}.    Bare mass slopes: {bareSlopes[0]:.2f}, {bareSlopes[1]:.2f}')

    with open('bareMassSlopes.in', 'w') as f:
        for bareSlope in bareSlopes:
            f.write(f'{bareSlope}\n')

    process = subprocess.call(f'make {programName} > /dev/null', shell=True)
    process = subprocess.call(f'./{programName} > /dev/null', shell=True)


    # ---------------------Plot the energy spectrum---------------------
    H_eigs_file = np.loadtxt(f"data/H_eigenvalues_{dataVar}_fit{fitNum}.out", skiprows=0)
    H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_{dataVar}_fit{fitNum}.out", skiprows=0)
    bare_data = np.loadtxt(f"data/bare_state_{dataVar}_fit{fitNum}.out", skiprows=0)
    bare_state_all = bare_data[:,0:4]
    bare_index_all = bare_data[:,4:7]
    info = np.loadtxt(f"data/finiteParams_{dataVar}_fit{fitNum}.out", skiprows=1)
    L = info[2]
    Lam_max = info[3]

    # Separate bare_index_all
    nr, nc = bare_state_all.shape
    bare_state = np.zeros((n_bare, int(nr/2), int(nc)))
    nr, nc = bare_index_all.shape
    bare_index = np.zeros((n_bare, int(nr/2), int(nc)))


    # CHECK THIS
    bare_state[1,:,:] = bare_state_all[0::2]
    bare_state[0,:,:] = bare_state_all[1::2]

    bare_index[1,:,:] = bare_index_all[0::2]
    bare_index[0,:,:] = bare_index_all[1::2]


    fig, ax = pypl.subplots()

    # ------------------Plot first odd parity nucleon-------------------
    ax.errorbar(m_pi_N2**2, m_N1, xerr=m_pi2_N1_err, yerr=m_N1_err,
                  fmt='.', ecolor='black', color='black', capsize=2
                  , label='CSSM', zorder=8.0)

    # ------------------Plot second odd parity nucleon------------------
    pypl.errorbar(m_pi_N2**2, m_N2, xerr=m_pi2_N2_err, yerr=m_N2_err,
                  fmt='.', ecolor='black', color='black', capsize=2
                  , label='_nolegend_', zorder=8.0)

    m_pi2 = H_eigs_file[:,0]
    n_m_pi = len(m_pi2)
    H_eigs = H_eigs_file[:,1:]
    nEigs = len(H_eigs[0,:])
    H0_eigs = H0_eigs_file[1:,1:]
    ch_num = H0_eigs_file[0,1:].astype(int)

    for i in range(0,nEigs):
        # Plot eigenvalues
        pypl.plot(m_pi2, H_eigs[:,i], color='black'
                  , label='_nolegend_', zorder=2.0)

    # Plot physical pion mass
    pypl.plot([m_pi0**2, m_pi0**2], [-4.0, 4.0], linestyle='dashed'
              , color='black', label='_nolegend')



    # ---------------------Bare state line plotting---------------------
    doPlotBare = True

    j = 0
    old_min1 = 0
    old_min2 = 0
    old_min3 = 0
    if doPlotBare:
        while (j < n_m_pi):
            if ((bare_index[0,j,0] != bare_index[0,j-1,0]) or (j==(n_m_pi-1))):
                pypl.plot(m_pi2[old_min1:(j)], bare_state[0,old_min1:(j),1]
                          , color='red', linewidth=3, alpha=1.0
                          , zorder=5.0, label='_nolegend')
                old_min1 = j
            if ((bare_index[0,j,1] != bare_index[0,j-1,1]) or (j==(n_m_pi-1))):
                pypl.plot(m_pi2[old_min2:(j)], bare_state[0,old_min2:(j),2]
                          , color='orange', linewidth=6, linestyle='dashed'
                          , alpha=1.0, zorder=3.0, label='_nolegend')
                old_min2 = j
            j = j + 1

    j = 0
    old_min1 = 0
    old_min2 = 0
    old_min3 = 0
    if doPlotBare:
        while (j < n_m_pi):
            if ((bare_index[1,j,0] != bare_index[1,j-1,0]) or (j==(n_m_pi-1))):
                pypl.plot(m_pi2[old_min1:(j)], bare_state[1,old_min1:(j),1]
                          , color='green', linewidth=3
                          , alpha=1.0, zorder=5.0, label='_nolegend')
                old_min1 = j
            if ((bare_index[1,j,1] != bare_index[1,j-1,1]) or (j==(n_m_pi-1))):
                pypl.plot(m_pi2[old_min2:(j)], bare_state[1,old_min2:(j),2]
                          , color='cyan', linewidth=6, linestyle='dashed'
                          , alpha=1.0, zorder=4.0, label='_nolegend')
                old_min2 = j
            j = j + 1

    # ------------------------------------------------------------------
    if doPlotBare:
        pypl.plot([-5.0,-5.0],[-5.0,-5.0], color='red'
                  , linewidth=3.0, label='Largest $m_{N_1}$')
        pypl.plot([-5.0,-5.0],[-5.0,-5.0], color='orange'
                  , linewidth=6.0, linestyle='dashed'
                  , label='2nd Largest $m_{N_1}$')
        pypl.plot([-5.0,-5.0],[-5.0,-5.0], color='green'
                  , linewidth=3.0, label='Largest $m_{N_2}$')
        pypl.plot([-5.0,-5.0],[-5.0,-5.0], color='cyan'
                  , linewidth=6.0, linestyle='dashed'
                  , label='2nd Largest $m_{N_2}$')
        pypl.legend(loc='lower right', numpoints=1)

    # ------------------------------------------------------------------

    if varL:
        pypl.title(f'S11 2b3c, Slopes: {bareSlopes[0]:.2f}, {bareSlopes[1]:.2f}')
    else:
        pypl.title(f'S11 2b3c, L = {L:.2f} fm, Slopes: {bareSlopes[0]:.2f}, {bareSlopes[1]:.2f}')
    max_E = 2.0
    pypl.axis([0.0, 0.4, 1.0, max_E])

    pypl.ylabel('E (GeV)')
    pypl.xlabel('$m_{\pi}^2 (\\textrm{GeV}^2)$')
    outFile = f'gridSearch/grid_2b3c_fit{fitNum}_{dataVar}_alpha_{int(bareSlopes[0]*100)}_{int(bareSlopes[1]*100)}.png'

    savefig(outFile, bbox_inches='tight')
    # pypl.show()

    pypl.close()
