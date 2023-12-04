#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math
import sys

# get fit number from file
with open('allFits.params', 'r') as f:
    f.readline()
    nFitLine = f.readline()
fitNum = int(nFitLine[-3:])

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

colours = ['red', 'blue']

# ------------------------------------------------------------------
plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, ax = plt.subplots()

# -------------------------bare state index-------------------------
bare_data = np.loadtxt(f"data/bare_state_lqcd_fit{fitNum}.out", skiprows=0)
bare_state_all = bare_data[:,0:4]
bare_index_all = bare_data[:,4:7]
# Separate bare_index_all
nr, nc = bare_state_all.shape
bare_state = np.zeros((n_bare, int(nr/2), int(nc)))
nr, nc = bare_index_all.shape
bare_index = np.zeros((n_bare, int(nr/2), int(nc)))

# CHECK THIS
bare_state[0,:,:] = bare_state_all[0::2]
bare_state[1,:,:] = bare_state_all[1::2]
bare_index[0,:,:] = bare_index_all[0::2]
bare_index[1,:,:] = bare_index_all[1::2]

# ------------------------------------------------------------------
cont_filename = f'data/contamination_lqcd_fit{fitNum}.out'
data_cont = np.loadtxt(cont_filename)
nMass = int(len(data_cont[:,0]) / 2)
nt = int(len(data_cont[0,:])) - 1

m_pi2 = np.zeros(nMass)
L = np.zeros(nMass)
t = np.zeros((nMass,nt))
Ct = np.zeros((nMass,nt))

for i in range(nMass):
    m_index = i*2
    L_index = i*2 + 1
    m_pi2[i] = data_cont[m_index,0]
    L[i] = data_cont[L_index,0]
    t[i,:] = data_cont[m_index,1:]
    Ct[i,:] = data_cont[L_index,1:]

Lfm = int(round(L[0]))

# choose which m_pi2 to plot
# 3 fm
# m_pi2_3fm = [
#     0.0286,
#     0.0784,
#     0.1529,
#     0.270,
#     0.388
# ]
m_pi2_3fm = [
    0.0290,
    0.0784,
    0.1610,
    0.270,
    0.388
]
# m_pi2_choice = 0.1385**2 # L&V point
# m_pi2_choice = 0.265 # wrong energy level
# m_pi2_choice = 0.2677 # also wrong now

# 2 fm
# m_pi2_2fm = [
#     0.266**2,
#     0.396**2,
#     0.444**2,
#     0.524**2
# ]
m_pi2_2fm = [
    0.15686,
    0.19711,
    0.27465
]
# m_pi2_choice = 0.185

massChoice = 4
if Lfm==2:
    m_pi2_choice = m_pi2_2fm[massChoice-1]
elif Lfm==3:
    m_pi2_choice = m_pi2_3fm[massChoice-1]


# m_pi2_choice = 0.2


doTitle = (True if (massChoice==0) else False)
# doTitle = (True if (massChoice==3) else False)

ind = np.argmin(np.abs(m_pi2 - m_pi2_choice))

redState = int(bare_index[0,ind,0])
blueState = int(bare_index[1,ind,0])
# print(f'red state: {redState}, blue state: {blueState}')

# ------------------------------------------------------------------
lineColour = 'red'
# lineColour = 'blue'

ab_choice = 'HEFT'
# ab_choice = 'lqcd'

eigsRemoved = 1
# eigsRemoved = 2
# ------------------------------------------------------------------

plotLabel =  '$m_{\pi}^{2}$' + f' = {m_pi2[ind]:.3f} GeV' + '${}^2$, ' \
    + '$E_{%i}$'  % ((redState if lineColour=='red' else blueState))

ax.plot(t[ind], Ct[ind], color=lineColour, label=plotLabel, linewidth=1.0)

yLabel = '$C_{%i}(t)$'  % ((1 if lineColour=='red' else 2))
ax.set_ylabel(yLabel, fontsize=26)
# ax.set_ylabel('$C(t)$')

if doTitle:
    # ax.set_title(f'S11 {int(n_bare)}b{int(n_ch)}c, ' \
        #              + f'L = {L[ind]:.1f} fm, ' \
        #              + '$m_{\pi}^{2}$' + f' = {m_pi2[ind]:.3f} GeV' + '${}^2$, ' \
        #              + '$E_{%i}$'  % ((redState if lineColour=='red' else blueState)))
    ax.set_title('$\\alpha_{j}$' + ',' + '$\\beta_{j}$ from ' \
                 + ('HEFT' if ab_choice=='HEFT' else 'lattice QCD'), fontsize=28)
ax.set_xlabel('$t$ (fm)', fontsize=26)

ax.set_xlim(0.0, L[ind])
ax.set_ylim(0.0, 0.8)

legendLoc = ('lower right' if Ct[ind,-1]>0.4 else 'upper right')

plt.legend(loc=legendLoc, numpoints=1
           , prop={'size': 24}, frameon=False)


# outfile = f'figs/S11_{int(n_bare)}b{int(n_ch)}c_contamination_varyingAB_' \
# outfile = f'figs/S11_{int(n_bare)}b{int(n_ch)}c_contamination_' \
#     + ab_choice + '_' + lineColour \
#     + f'_L{int(L[ind]*10)}_mpisq{int(m_pi2[ind]*1000)}_fit{fitNum}.{figType}'

# outfile = f'figs/S11_{int(n_bare)}b{int(n_ch)}c_contamination_' \
#     + ab_choice + '_' + lineColour \
#     + f'_L{int(L[ind]*10)}_mpisq{int(m_pi2[ind]*1000)}.{figType}'
outfile = f'figs/S11_{int(n_bare)}b{int(n_ch)}c_contamination_' \
    + ab_choice + '_' + lineColour + f'_{int(eigsRemoved)}eigs'\
    + f'_L{Lfm}fm_M{int(massChoice)}.{figType}'


# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')

if saveFigs:
    savefig(outfile, bbox_inches='tight')
else:
    plt.show()
