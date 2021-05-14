##################################################
# This script plots phase diagram of a CW model. #
##################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10., 6.))

plt.ylim([0,3])
plt.xlim([-1,1])

beta = 1
J = np.arange(1,3,.0001)
h1 = -J*np.sqrt(1-1/(beta*J)) + 1/(2*beta)*( np.log(1+np.sqrt(1-1/(beta*J))) - np.log(1-np.sqrt(1-1/(beta*J))) )
h2 = J*np.sqrt(1-1/(beta*J)) + 1/(2*beta)*( np.log(1-np.sqrt(1-1/(beta*J))) - np.log(1+np.sqrt(1-1/(beta*J))) )

plt.xlabel(r'$\beta h$', fontsize=20)
plt.ylabel(r'$\beta J$', fontsize=20)

plt.axvline(0, alpha=.5, color='black')
plt.axhline(0, alpha=.5, color='black')

plt.plot(h1, J, color='red', label='Singular curve')
plt.plot(h2, J, color='red')
plt.plot(0,1, marker="+", markersize=15, color='blue', linestyle='none', mew=1.5, alpha=.9, label="Critical point")

plt.fill_between(h1, J, 1000, color='springgreen', alpha=.3, label='Bistable region')
plt.fill_between(h2, J, 1000, color='springgreen', alpha=.3)

plt.grid()
plt.legend(loc='lower right', prop={'size': 20})

############### Zoomed insert #################
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

axins = zoomed_inset_axes(ax, 10, bbox_to_anchor=(0.05, 0.05, .6, .5), bbox_transform=ax.transAxes, loc="lower left")
axins.axhline(0, color='black', alpha=.5)
axins.axvline(0, color='black', alpha=.5)

axins.plot(h1, J, color='red', label='Singular curve')
axins.plot(h2, J, color='red')
axins.plot(0,1, marker='+', markersize=15, mew=1.5, linestyle='none', color='blue', zorder=10, alpha=.9)

axins.fill_between(h1, J, 1000, color='springgreen', alpha=.3, label='Bistable region')
axins.fill_between(h2, J, 1000, color='springgreen', alpha=.3)

x1, x2, y1, y2 = -.03, .03, .99, 1.09 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
plt.grid()

mark_inset(ax, axins, loc1=1, loc2=4, fc="none", ec="0.5", zorder=0)

plt.tight_layout()
plt.savefig("../figures/CW_phase_diagram.png", dpi=300)
plt.show()
