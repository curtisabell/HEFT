#!/usr/bin/python3
import numpy as np
from pylab import *
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
import math

plt.style.use('classic')
plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
fig, ax = plt.subplots()

# nucleon stuff
# a_fm = [0.1687]


m_hadron = np.array([1.0132855, 1.0638650, 1.1525865, 1.2651921, 1.4091008])
dm_hadron = np.array([0.0049499, 0.0053206, 0.0058410, 0.0127765, 0.0371326])

m_pi0 = 0.1385
m_hadron2 = m_hadron**2
m_pi = np.array([0.1687, 0.2799, 0.3905, 0.5149, 0.6228])
m_pi2 = m_pi**2
dm_pi = m_pi - m_pi0
dm_pi2 = m_pi**2 - m_pi0**2


theta_N = np.polyfit(dm_pi2, m_hadron, 1)
m_pi_fit = np.array([0.1385, 0.1687, 0.2799, 0.3905, 0.5149, 0.6228])
m_pi2_fit = m_pi_fit**2
dm_pi2_fit = m_pi2_fit - m_pi0**2
m_hadron_fit = theta_N[1] + theta_N[0]*dm_pi2_fit


# Plot in m_pi
ax.errorbar(m_pi, m_hadron, yerr=dm_hadron,
            fmt='.', ecolor='black', color='black'
            , capsize=2, label='_nolegend_', zorder=10.0)
# ax.scatter(m_pi, m_hadron, color='red')
ax.plot(m_pi_fit, m_hadron_fit, color='red')
# ax.plot(m_pi, m_hadron, color='black', linestyle='dashed')

# # Plot in m_pi2
# ax.scatter(m_pi2, m_hadron, color='red')
# ax.plot(m_pi2_fit, m_hadron_fit)
# ax.plot(m_pi2, m_hadron, color='red', linestyle='dashed')




# ax.set_title(f'mN0 = {theta_N[1]:.4f}, slope = {theta_N[0]:.4f}')
ax.set_title('$m_{N}$ = '+f'{theta_N[1]:.4f} + {theta_N[0]:.4f}' + '$(m_{\pi}^2 - m_{\pi}^2|_{phys})$')

# theta_N = np.polyfit(dm_pi2, m_hadron2, 1)
# dm_pi2_fit = m_pi2_fit - m_pi0**2


# ax.scatter(m_pi2, m_hadron)
# E_piN_0 = m_pi + m_hadron
# # ax.plot(m_pi2, E_piN_0, color='red', linestyle='solid')
# ax.scatter(m_pi2, E_piN_0, color='red')
# ax.set_ylim([1.0, 2.05])
# ax.set_xlim([0.0, 0.4])

# E_piN_0_fit1 = np.sqrt(m_pi2_fit) + sqrt(m_hadron2_fit)
# ax.plot(m_pi2_fit, E_piN_0_fit1, color='red', linestyle='dashed')
# # E_piN_0_fit1 = np.sqrt(m_pi2_fit) + np.sqrt(theta_N[1]) + theta_N[0]*m_pi2_fit



# print(thetaConv(theta_N[1], theta_N[0]))

# xtickloc = []
# ax.xticks(xtickloc)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', width=1, length=10, color='black')
ax.tick_params(axis='x', which='major', width=1, length=10, color='black')
ax.tick_params(axis='y', which='minor', width=1, length=4,  color='black', direction='in', bottom='on')
ax.tick_params(axis='x', which='minor', width=1, length=4,  color='black', direction='in', top='on')
plt.show()
