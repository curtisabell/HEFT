#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as pypl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

sys.path.append('../src')
import readHEFTConfig
HEFT = readHEFTConfig.HEFTConfig()
HEFT.readHEFTConfigFile()

pypl.style.use('classic')
pypl.rcParams['axes.formatter.useoffset'] = False
pypl.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
pypl.rc('text', usetex=True)
fig, ax = pypl.subplots()

m_pi0 = 0.1385
# ---------------------------Get fit info---------------------------
# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

# ------------------------Check to save figs------------------------
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

# --------------------Read mpiFiniteVol.f90 data--------------------
with open(f'data/latest_finiteProgram_fit{fitNum}.out', 'r') as f:
    finProg = f.readline().strip()
if (finProg=='lqcd'):
    plotVarL = True
    programString = 'lqcd'
elif (finProg=='mpi'):
    plotVarL = False
    programString = 'm_pi'
else:
    print('Error in latest finiteProgram file')
    exit()

H_eigs_file = np.loadtxt(f"data/H_eigenvalues_{programString}_fit{fitNum}.out", skiprows=0)
H0_eigs_file = np.loadtxt(f"data/H0_eigenvalues_{programString}_fit{fitNum}.out", skiprows=0)
bare_data = np.loadtxt(f"data/bare_state_{programString}_fit{fitNum}.out", skiprows=0)
info = np.loadtxt(f"data/finiteParams_{programString}_fit{fitNum}.out", skiprows=1)

bare_state_all = bare_data[:,0:4]
bare_index_all = bare_data[:,4:7]
n_ch = int(info[0])
n_bare = int(info[1])
L = info[2]
Lfm = int(round(L))
Lam_max = info[3]
print(Lfm)

m_pi_D200 = 0.2
lqcd_capsize=4.0
lqcd_linewidth=2.0
lqcd_fmt = 'o'
doPlotLQCD = True
if doPlotLQCD:
    # Morningstar 200 MeV results
    #G_1u(0):
    ax.errorbar([m_pi_D200], [1.15], yerr=[0.009],
                fmt=lqcd_fmt, ecolor='black', color='black', capsize=lqcd_capsize
                , label='D200', zorder=8.0, elinewidth=lqcd_linewidth)
    ax.errorbar([m_pi_D200], [1.3565], yerr=[0.0125],
                fmt=lqcd_fmt, ecolor='black', color='black', capsize=lqcd_capsize
                , label='_nolegend_', zorder=8.0, elinewidth=lqcd_linewidth)
    # G_1(1):
    # ax.errorbar([m_pi_D200], [1.158], yerr=[0.009],
    #             fmt='.', ecolor='black', color='black', capsize=2
    #             , label='D200', zorder=8.0)
    # ax.errorbar([m_pi_D200], [1.279], yerr=[0.010],
    #             fmt='.', ecolor='black', color='black', capsize=2
    #             , label='_nolegend_', zorder=8.0)
# ------------------------------------------------------------------




# -------------------------Plot mpiFin data-------------------------
m_pi2 = H_eigs_file[:,0]
n_m_pi = len(m_pi2)
H_eigs = H_eigs_file[:,1:]
nEigs = len(H_eigs[0,:])
H0_eigs = H0_eigs_file[1:,1:]
ch_num = H0_eigs_file[0,1:].astype(int)

longdash = (0, (12, 8))
longdashdot = (0, (12, 8, 2, 8))

ch_labels = ['Bare mass', '$\pi N$', '$\eta N$', '$K\Lambda$']
# basisColors = ['gray', 'blue', 'red', 'green']
basisColors = ['gray', 'orange', 'green','darkviolet']
style1 = 'dashed'
style2 = (0, (14, 8))
style3 = (0, (8, 6, 3, 6))
basisStyles = ['solid', style2, style1, style3]
ch_style = [longdashdot, 'dashed', 'dashed', 'dashed']
lwidth = 1.0

plotBasis = True
plotEigs = True
plotBareBasis = False

D200_index = np.argmin(np.abs(m_pi_D200**2 - m_pi2))

# Plot bare state first
if plotBasis and plotBareBasis:
    pypl.plot(np.sqrt(m_pi2), H0_eigs[:,0], color='gray'
              , linestyle=longdashdot, label='_nolegend_'
              , zorder=1.0, linewidth=lwidth)
    pypl.plot(np.sqrt(m_pi2), H0_eigs[:,1], color='gray'
              , linestyle=longdashdot, label='_nolegend_'
              , zorder=1.0, linewidth=lwidth)


for i in range(0,nEigs):
    # Plot free-particle states
    if plotBasis:
        if i>1:
            pypl.plot(np.array([0.0,1.0]), H0_eigs[D200_index,i]*np.array([1,1])
                      , color=basisColors[ch_num[i]], zorder=2.0
                      , linewidth=lwidth
                      , linestyle=basisStyles[ch_num[i]], label='_nolegend_')

    if plotEigs:
        # Plot eigenvalues
        pypl.plot(np.array([0.0,1.0]), H_eigs[D200_index,i]*np.array([1,1])
                  , color='black'
                  , label='_nolegend_', zorder=2.0)


# ---------------------Bare state contributions---------------------
doPlotBare = False
if doPlotBare:
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

    bs_colours = ['blue', 'red']

    j = 0
    old_min1 = 0
    old_min2 = 0
    old_min3 = 0
    if doPlotBare:
        while (j < n_m_pi):
            if ((bare_index[0,j,0] != bare_index[0,j-1,0]) or (j==(n_m_pi-1))):
                pypl.plot(np.sqrt(m_pi2)[old_min1:(j)], bare_state[0,old_min1:(j),1]
                          , color=bs_colours[0], linewidth=3, alpha=1.0
                          , zorder=5.0, label='_nolegend')
                old_min1 = j
            if ((bare_index[0,j,1] != bare_index[0,j-1,1]) or (j==(n_m_pi-1))):
                pypl.plot(np.sqrt(m_pi2)[old_min2:(j)], bare_state[0,old_min2:(j),2]
                          , color=bs_colours[0], linewidth=4, linestyle=longdash
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
                pypl.plot(np.sqrt(m_pi2)[old_min1:(j)], bare_state[1,old_min1:(j),1]
                          , color=bs_colours[1], linewidth=3
                          , alpha=1.0, zorder=4.0, label='_nolegend')
                old_min1 = j
            if ((bare_index[1,j,1] != bare_index[1,j-1,1]) or (j==(n_m_pi-1))):
                pypl.plot(np.sqrt(m_pi2)[old_min2:(j)], bare_state[1,old_min2:(j),2]
                          , color=bs_colours[1], linewidth=4, linestyle=longdash
                          , alpha=1.0, zorder=4.0, label='_nolegend')
                old_min2 = j
            j = j + 1


# ---------------------------Legend stuff---------------------------
plotBareLegend = False
if doPlotBare and plotBareLegend:
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=bs_colours[0]
              , linewidth=3.0, label='Largest $m_{N_1}$')
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=bs_colours[0]
              , linewidth=4.0, linestyle=longdash
              , label='2nd Largest $m_{N_1}$')
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=bs_colours[1]
              , linewidth=3.0, label='Largest $m_{N_2}$')
    pypl.plot([-5.0,-5.0],[-5.0,-5.0], color=bs_colours[1]
              , linewidth=4.0, linestyle=longdash
              , label='2nd Largest $m_{N_2}$')
    pypl.legend(loc='lower right', numpoints=1
                , prop={'size': 18}, frameon=False)


# plotBasisLegend = True
# if plotBasis and plotBasisLegend and not (doPlotBare and plotBareLegend):
#     for i in range(4):
#         if (i==0 and (not plotBareBasis)):
#             continue
#         pypl.plot([-5.0,-5.0],[-5.0,-5.0]
#                   , color=basisColors[i], linestyle=basisStyles[i]
#                   , linewidth=lwidth, label=ch_labels[i])
#     pypl.legend(loc='lower right', numpoints=1
#                 , prop={'size': 14}, frameon=False)


pypl.plot([-5.0,-5.0],[-5.0,-5.0]
          , color=basisColors[1], linestyle=basisStyles[1]
          , linewidth=lwidth, label=ch_labels[1])
pypl.legend(loc='lower right', numpoints=1
            , prop={'size': 22}, frameon=False)



# Plot physical pion mass
pypl.plot([m_pi0**2, m_pi0**2], [-4.0, 4.0], linestyle='dashed'
          , color='black', label='_nolegend')

if plotVarL:
    # outFile = f'figs/2b3c_Ev{programString}_fit{fitNum}.{figType}'
    # titleString = f'S11 N2b3c, Fit {fitNum}'
    titleString = 'L ' + '$\sim$' + f' {int(L)} fm'
    outFile = f'figs/2b3c_EvMpi_basisOnly_{int(L)}fm.{figType}'
else:
    outFile = f'figs/2b3c_Ev{programString}_L_{int(L*10)}_fit{fitNum}.{figType}'
    titleString = f'S11 N2b3c, L = {L:.2f} fm,  Fit {fitNum}'

outFile = f'figs/2b3c_EvMpi_D200.{figType}'
titleString = '$L$' + f' = {L:.2f} fm'
# pypl.title(titleString)

max_E = 1.4
x_fac = m_pi_D200*0.01
# pypl.axis([0.0, 0.4, 1.0, max_E])
# pypl.axis([0.0, max(np.sqrt(m_pi2)), 1.0, max_E])
# pypl.axis([0.038, 0.042, 1.0, max_E])
pypl.axis([m_pi_D200 - x_fac, m_pi_D200 + x_fac, 1.0, max_E])

pypl.ylabel('E (GeV)', fontsize=22)
pypl.xlabel('$m_{\pi} (\\textrm{GeV})$', fontsize=22)

# xtickloc = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
# pypl.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
pypl.tick_params(axis='y', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='x', which='major', width=1, length=10, color='black')
pypl.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
pypl.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

fig.set_size_inches(5, 8)
pypl.xticks([m_pi_D200])

if not saveFigs:
    pypl.show()
else:
    savefig(outFile, bbox_inches='tight')
