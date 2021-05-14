##############################################
# This script produces plots that illustrate #
# graphical analysis for CW model.           #
##############################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(23/2/1.1, 10.5/1.2))

T = 1

def A(m):
    return T/2*np.log((1+m)/(1-m))

x_artanh = np.arange(-1,1,.00001) #.00001
y_artanh = A(x_artanh)

J, h = .5, 0
title = r"$\beta J=$" + str(J) + r"$,\ \beta h=$" + str(h)
axs[0,0].set_xlim([-1.015, 1.015])
axs[0,0].set_ylim([-2, 2])
axs[0,0].set_ylabel(r'${\cal A}(\bar m),\ {\cal P}(\bar m; J,h,\beta)$', fontsize=18)
axs[0,0].axhline(0, color='black')
axs[0,0].axvline(0, color='black')
X = np.arange(-1, 1, .001)

def P(m):
    return J*m + h
def SC(m):
    return P(m)-A(m)

m_stab_ = fsolve(SC, 0)[0]
axs[0,0].plot(X, P(X), color='green', zorder=11, linewidth=2.5, label=r'${\cal P}(\bar m; J,h,\beta)$')
axs[0,0].plot(x_artanh, y_artanh, color='blue', zorder=10, linewidth=2.5, label=r'${\cal A}(\bar m)$')
axs[0,0].plot(m_stab_, P(m_stab_), marker='^', linestyle='none', color='red', markersize=7, label='Stable solution', zorder=100)

axs[0,0].legend(loc='lower right', fontsize=11.5)
axs[0,0].set_title(title, fontsize=18)
axs[0,0].grid()

####################################################################################
axs[0,1].set_xlim([-1.015, 1.015])
axs[0,1].set_ylim([-2, 2])
J, h = 2, 0
title = r"$\beta J=$" + str(J) + r"$,\ \beta h=$" + str(h)
axs[0,1].axhline(0, color='black')
axs[0,1].axvline(0, color='black')

m_stab = np.array([fsolve(SC, -.99)[0], fsolve(SC, .99)[0]])
m_unstab = fsolve(SC, 0)[0]
axs[0,1].plot(X, P(X), color='green', zorder=11, linewidth=2.5, label=r'${\cal P}(\bar m; J,h,\beta)$')
axs[0,1].plot(x_artanh, y_artanh, color='blue', zorder=10, linewidth=2.5, label=r'${\cal A}(\bar m)$')
axs[0,1].plot(m_stab, P(m_stab), marker='^', linestyle='none', color='red', markersize=7, label='Stable solutions', zorder=100)
axs[0,1].plot(m_unstab, P(m_unstab), marker='v', linestyle='none', color='yellow', markersize=7, label='Unstable solution', zorder=100)

axs[0,1].legend(loc='lower right', fontsize=11.5)

axs[0,1].set_title(title, fontsize=18)
axs[0,1].grid()

##################### Free energy ###########################
J, h = .5, 0
dm = .001
m = np.arange(-1, 1+dm, dm)
def F(m_):
    return -J*m_**2/2 - h*m_ + T/2*((1+m_)*np.log((1+m_)) + (1-m_)*np.log(1-m_)) - np.log(2)
axs[1,0].plot(m, F(m), linewidth=3, label='SMFE')
axs[1,0].plot(m_stab_, F(m_stab_), marker='^', linestyle='none', color='red', markersize=7, label='Stable solution', zorder=100)
axs[1,0].set_xlim([-1.015, 1.015])
axs[1,0].set_ylim([-.7, -.24])
axs[1,0].set_xlabel(r'$\bar m$', fontsize=18)
axs[1,0].set_ylabel(r'$\bar{\cal F}(\bar m; J,h,\beta)$', fontsize=18)
axs[1,0].grid()
axs[1,0].axvline(0, color='black')
axs[1,0].axhline(0, color='black')
axs[1,0].legend(loc=[.54, .83], fontsize=11.5)

#############################################################
J, h = 2, 0
axs[1,1].plot(m, F(m), linewidth=3, label='SMFE')
axs[1,1].plot(m_stab, F(m_stab), marker='^', linestyle='none', color='red', markersize=7, label='Stable solutions', zorder=100)
axs[1,1].plot(m_unstab, F(m_unstab), marker='v', linestyle='none', color='yellow', markersize=7, label='Unstable solution', zorder=100)
axs[1,1].set_xlim([-1.015, 1.015])
axs[1,1].set_ylim([-1.05, -.68])
axs[1,1].set_xlabel(r'$\bar m$', fontsize=18)
axs[1,1].grid()
axs[1,1].axvline(0, color='black')
axs[1,1].axhline(0, color='black')
axs[1,1].legend(loc='lower center', fontsize=11.5)

plt.tight_layout()
plt.savefig('../figures/CW_graphical_analysis_with_FE.png', dpi=300)
plt.show()

fig.clf()
plt.close()
