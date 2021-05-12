#############################################################
# This script plots the specific microcanonical free energy #
# of the MCW model in the vicinity of a critical point.     #
#############################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

gamma1 = .5
gamma2 = 1-gamma1

J11 = 1
J22 = 0.487664
J12 = -1.615554

h2 = 0.669723
h1 = 0.563408


def F(m):
    return -.5*gamma1**2*J11*m[0]**2 - .5*gamma1**2*J22*m[1]**2 - gamma1*gamma2*J12*m[0]*m[1] - gamma1*h1*m[0] - gamma2*h2*m[1] + gamma1/2*((1-m[0])*np.log(1-m[0])+(1+m[0])*np.log(1+m[0])) - gamma1*np.log(2) + gamma2/2*((1-m[1])*np.log(1-m[1])+(1+m[1])*np.log(1+m[1])) - gamma2*np.log(2)

X = np.arange(.3-.45, .3+.4, 0.001)
xlen = len(X)
Y = np.arange(.5-.45, .5+.4, 0.001)
ylen = len(Y)
X, Y = np.meshgrid(X, Y)
Z = F([X, Y])

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False, vmin=-.822, vmax=-.821, alpha=.5)

ax.scatter(.3, .5, F([.3,.5]), color='darkviolet', marker='+', s=500, edgecolor='darkviolet', linewidth=3)

# Customize the z axis
ax.set_zlim(-.823, -.8)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.04f'))

ax.set_xlabel(r'$\bar m_1$', fontsize=35, labelpad=20)
ax.set_ylabel(r'$\bar m_2$', fontsize=35, labelpad=20)
ax.set_zlabel(r'$\bar{\cal F}$', fontsize=35, labelpad=22)

ax.tick_params(axis='x', labelsize=20)   
ax.tick_params(axis='y', labelsize=20)  
ax.tick_params(axis='z', labelsize=20)

plt.tight_layout()

plt.show()

