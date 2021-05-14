##########################################
# This script plots family of CW model's #
# specific microcanonical free energies. #
##########################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

from scipy.optimize import fsolve
import matplotlib.pyplot as plt 
import math
import numpy as np

fig, axs = plt.subplots(1,2, figsize=(18.5/1.2, 10.5/2))

J_ = [0.9, 0.95, 1, 1.05, 1.1]
T = 1.
colors = ['orange', 'r', 'g', 'cyan', 'b']

##### Plotting free energy ######
axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')

i = 0
m = np.arange(-1, 1, 0.0001)
h = 0
T = 1
for J in J_:
    F = -.5*J*m**2 - h*m + T/2*((1+m)*np.log((1+m)) + (1-m)*np.log(1-m))
    axs[0].plot(m, F, label=r'$\beta J=$' + str(J), linewidth=3, color=colors[i])
    i += 1
axs[0].set_xlim([-1, 1])
axs[0].set_ylim([-.012, .01])
axs[0].grid()
axs[0].set_xlabel(r'$\bar m$', fontsize=20)
axs[0].set_ylabel(r'$\bar{\cal F}(\bar m;J,h,\beta)$', fontsize=20, labelpad=-1)
axs[0].set_title(r'$\beta h=0;\ \beta J=\{0.9, 0.95, 1, 1.05, 1.1\}$', fontsize=18)
axs[0].legend(loc=[.37,.01], prop={'size': 13.1})

##### Plotting free energy ######
axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')

i = 0
m = np.arange(-1, 1, 0.0001)
h = .01
T = 1
for J in J_:
    F = -.5*J*m**2 - h*m + T/2*((1+m)*np.log((1+m)) + (1-m)*np.log(1-m))
    axs[1].plot(m, F, label=r'$\beta J=$' + str(J), linewidth=3, color=colors[i])
    i += 1
axs[1].set_xlim([-1, 1])
axs[1].set_ylim([-.012, .01])
axs[1].grid()
axs[1].set_xlabel(r'$\bar m$', fontsize=20)
axs[1].set_title(r'$\beta h=0.01;\ \beta J=\{0.9, 0.95, 1, 1.05, 1.1\}$', fontsize=18)
axs[1].legend(loc=[.37,.01], prop={'size': 13.1})

plt.tight_layout()
plt.savefig('../figures/CW_free_energies.png', dpi=300)
plt.show()

fig.clf()
plt.close()
