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

J, h = 2, -0.53283998
title = r"$\beta J=$" + str(J) + r"$,\ \beta h=$" + str(round(h,5))
axs[0,1].set_xlim([-1.015, 1.015])
axs[0,1].set_ylim([-2.7, 1.5])
axs[0,1].axhline(0, color='black')
axs[0,1].axvline(0, color='black')
X = np.arange(-1, 1, .001)

def P(m):
    return J*m + h
def SC(L):
    return P(L)-A(L)

m_stab = fsolve(SC, -.99)[0]
m_sing = fsolve(SC, .99)[0]
axs[0,1].plot(x_artanh, y_artanh, color='blue', zorder=10, linewidth=2.5, label=r'${\cal A}(\bar m)$')
axs[0,1].plot(X, P(X), color='green', zorder=11, linewidth=2.5, label=r'${\cal P}(\bar m; J, h,\beta)$')
axs[0,1].plot(m_sing, P(m_sing), marker='D', linestyle='none', color='orange', markersize=7, label='Spinodal point', zorder=100)
axs[0,1].plot(m_stab, P(m_stab), marker='^', linestyle='none', color='red', markersize=7, label='Stable solution', zorder=100)

axs[0,1].legend(loc='lower right', fontsize=12.5)
axs[0,1].set_title(title, fontsize=18)
axs[0,1].grid()
####################################################################################
J, h = 1, 0
title = r"$\beta J=$" + str(J) + r"$,\ \beta h=$" + str(h)
axs[0,0].set_xlim([-1.015, 1.015])
axs[0,0].set_ylim([-2.7, 1.5])
axs[0,0].set_ylabel(r'${\cal A}(\bar m),\ {\cal P}(\bar m; J, h,\beta)$', fontsize=18)
axs[0,0].axhline(0, color='black')
axs[0,0].axvline(0, color='black')

m_crit = fsolve(SC, 0)[0]
axs[0,0].plot(x_artanh, y_artanh, color='blue', zorder=10, linewidth=2.5, label=r'${\cal A}(\bar m)$')
axs[0,0].plot(X, P(X), color='green', zorder=11, linewidth=2.5, label=r'${\cal P}(\bar m; J, h,\beta)$')
axs[0,0].plot(m_crit, P(m_crit), marker='*', linestyle='none', color='red', markersize=9, label='Critical point', zorder=100)

axs[0,0].legend(loc='lower right', fontsize=12.5)
axs[0,0].set_title(title, fontsize=18)
axs[0,0] .grid()

##################### Free energy ###########################
J, h = 1, 0
dm = .001
m = np.arange(-1, 1+dm, dm)
def F(m_):
    return -J*m_**2/2 - h*m_ + T/2*((1+m_)*np.log((1+m_)) + (1-m_)*np.log(1-m_)) - np.log(2)
axs[1,0].plot(m, F(m), linewidth=3, label='SMFE')
axs[1,0].plot(m_crit, F(m_crit), marker='*', linestyle='none', color='red', markersize=9, label='Critical point', zorder=100)
axs[1,0].set_xlim([-1.015, 1.015])
axs[1,0].set_ylim([-.7, -.5])
axs[1,0].set_yticks(np.arange(-0.7, -.5+0.1, 0.1))
axs[1,0].set_xlabel(r'$\bar m$', fontsize=18)
axs[1,0].set_ylabel(r'$\bar{\cal F}(\bar m; J,h,\beta)$', fontsize=18)
axs[1,0].grid()
axs[1,0].axvline(0, color='black')
axs[1,0].axhline(0, color='black')
axs[1,0].legend(loc=[.545, .82], fontsize=12.5)

#############################################################
J, h = 2, -0.53283998
axs[1,1].plot(m, F(m), linewidth=3, label='SMFE')
axs[1,1].plot(m_sing, F(m_sing), marker='D', linestyle='none', color='orange', markersize=7, label='Spinodal point', zorder=100)
axs[1,1].plot(m_stab, F(m_stab), marker='^', linestyle='none', color='red', markersize=7, label='Stable solution', zorder=100)
axs[1,1].set_xlim([-1.015, 1.015])
axs[1,1].set_ylim([-1.6, -.4])
axs[1,1].set_xlabel(r'$\bar m$', fontsize=18)
axs[1,1].grid()
axs[1,1].axvline(0, color='black')
axs[1,1].axhline(0, color='black')
axs[1,1].legend(loc='lower right', fontsize=12.5)

plt.tight_layout()
plt.savefig('../figures/CW_graphical_analysis_singular_with_FE.png', dpi=300)
plt.show()

fig.clf()
plt.close()
