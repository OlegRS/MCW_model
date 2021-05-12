####################################################
# This script shows the plots needed for graphical #
# analysis of the 2-component CW model.            #
####################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(23/2, 10.5/2))
gamma1 = .5
gamma2 = .5
J11 = 3
J22 = 3
J12 = -.42
h1 = 0
h2 = 0

M1 = np.arange(-1, 1,.00001)
M2 = np.arange(-1, 1,.00001)

m2 = 1/(gamma2*J12)*(1/2*np.log((1+M1)/(1-M1)) - J11*gamma1*M1 - h1)
m1 = 1/(gamma1*J12)*(1/2*np.log((1+M2)/(1-M2)) - J22*gamma2*M2 - h2)

axs[0].set_xlabel(r'$\bar m_1$', fontsize=15)
axs[0].set_ylabel(r'$\bar m_2$', fontsize=15, labelpad=-5)
axs[0].axis('equal')
axs[0].set_xlim([-1, 1])
axs[0].set_ylim([-1, 1])
axs[0].axhline(0, linewidth=1, color='black')
axs[0].axvline(0, linewidth=1, color='black')
axs[0].plot(m2, M1, label=r'Equation_1: $\frac{\partial F}{\partial m_1} = 0$', linewidth=2)
axs[0].plot(M2, m1, label=r'Equation_2: $\frac{\partial F}{\partial m_2} = 0$', linewidth=2)

axs[0].set_title("$\gamma_1=$" + str(gamma1) + ", $\gamma_2=$" + str(gamma2) + "; $J_{11}=$"+str(J11)+", $J_{22}=$"+str(J22)+", $J_{12}=$"+str(J12)+"; $h_1=$"+str(h1)+", $h_2=$"+str(h2), fontsize=12, y=1.01)
axs[0].grid()

#########################################################################
gamma1 = .5
gamma2 = .5
J11 = 3
J22 = 3
J12 = -.42
h1 = 0
h2 = .2

M1 = np.arange(-1, 1,.00001)
M2 = np.arange(-1, 1,.00001)

m2 = 1/(gamma2*J12)*(1/2*np.log((1+M1)/(1-M1)) - J11*gamma1*M1 - h1)
m1 = 1/(gamma1*J12)*(1/2*np.log((1+M2)/(1-M2)) - J22*gamma2*M2 - h2)

axs[1].set_xlabel(r'$\bar m_1$', fontsize=15)
axs[1].axis('equal')
axs[1].set_xlim([-1, 1])
axs[1].set_ylim([-1, 1])
axs[1].axhline(0, linewidth=1, color='black')
axs[1].axvline(0, linewidth=1, color='black')
axs[1].plot(m2, M1, label=r'Equation_1: $\frac{\partial F}{\partial m_1} = 0$', linewidth=2)
axs[1].plot(M2, m1, label=r'Equation_2: $\frac{\partial F}{\partial m_2} = 0$', linewidth=2)

axs[1].set_title("$\gamma_1=$" + str(gamma1) + ", $\gamma_2=$" + str(gamma2) + "; $J_{11}=$"+str(J11)+", $J_{22}=$"+str(J22)+", $J_{12}=$"+str(J12)+"; $h_1=$"+str(h1)+", $h_2=$"+str(h2), fontsize=12, y=1.01)
axs[1].grid()

plt.tight_layout()
plt.savefig("../figures/2CW_graphical_analysis.png", dpi=300)
plt.show()
fig.clf()
plt.close()
