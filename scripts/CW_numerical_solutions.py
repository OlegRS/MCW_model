############################################
# This script plots family of solutions of #
# of the CW model's equation of state.     #
############################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

from scipy.optimize import fsolve
import matplotlib.pyplot as plt 
import math
import numpy as np

J_ = [0, .5, 1, 1.5, 2]
T = 1.
colors = ['orange', 'r', 'g', 'cyan', 'b']

##### Plotting self-consistency equation ######
fig = plt.figure(dpi=200)
plt.axvline(0, color='black', alpha=.8)
plt.axhline(0, color='black', alpha=.8)

i = 0
for J in J_:
    delta = 0.0025
    Frange = np.arange(-1.5, 1.5, delta)
    mrange = np.arange(-1.1, 1.1, delta)
    M, F = np.meshgrid(mrange, Frange)
    G = M - np.tanh(J/T*M + F/T)
    CS = plt.contour(F, M, G, [0], colors=colors[i], linewidths=1.5)
    CS.collections[0].set_label(r'$\beta J=$' + str(J))
    i += 1
plt.grid()
plt.xlabel(r'$\beta h$', fontsize=15)
plt.ylabel(r'$\bar m$', fontsize=15, labelpad=-3)
plt.ylim([-1.01, 1.01])
plt.legend(loc=4, prop={'size': 10})
fig.set_size_inches(18.5/3, 10.5/3)

plt.tight_layout()
plt.savefig('../figures/CW_numerical_solutions.png', dpi=300)
plt.show()

fig.clf()
plt.close()
