####################################################
# This script plots magnetisations as functions of #
# field indicating stable solutions.               #
####################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
import random

gamma1 = .5
gamma2 = 1-gamma1

J11 = 3
J22 = 3
J12 = -.42

h2_l = -.4
h2_r = .4
h2_step = .001

title = r"$J_{11}$="+str(J11)+ r"$,\ J_{22}=$"+str(J22)+r"$,\ J_{12}=$"+str(J12) + r"$;\ \gamma_1=$" + str(gamma1) + r"$,\ \gamma_2=$" + str(gamma2)#+ r"$;\ h_1=$"+str(h1)

def equations_of_state(m):
    eq1 = h1 + J11*gamma1*m[0] + J12*gamma2*m[1] - 1/2*np.log((1+m[0])/(1-m[0]))
    eq2 = h2 + J22*gamma2*m[1] + J12*gamma1*m[0] - 1/2*np.log((1+m[1])/(1-m[1]))
    return (eq1, eq2)

def eigenvalues(m):
    J11_ = J11-1/(gamma1*(1-m[0]**2))
    J22_ = J11-1/(gamma2*(1-m[1]**2))
    lambda1 = -(J11_ + J22_ + np.sqrt((J11_-J22_)**2 + 4*J12**2))/2
    lambda2 = -(J11_ + J22_ - np.sqrt((J11_-J22_)**2 + 4*J12**2))/2
    return lambda1, lambda2

def Jacobian_matrix(m):
    return ((J11*gamma1-1/(1-m[0]**2), J12*gamma2), (J12*gamma1, J22*gamma2-1/(1-m[1]**2)))

# Scanning h2
H2 = []
M1 = []
M2 = []
H2_stable = []
M1_stable = []
M2_stable = []

bnds = ((-.99999, .99999), (-.99999, .99999))
for h2 in np.arange(h2_l, h2_r+h2_step, h2_step):
    h1 = -h2
    completeness = (h2-h2_l)/(h2_r-h2_l)*100
    if completeness % 10 <= h2_step/(h2_r-h2_l)*100:
        print(int(completeness), "% complete")
    for count in range(0, 300):
        m1, m2 = optimize.fsolve(equations_of_state, [random.uniform(bnds[0][0], bnds[0][1]), random.uniform(bnds[1][0], bnds[1][1])], fprime = Jacobian_matrix)
        eq1, eq2 = equations_of_state([m1, m2])
        if eq1**2 + eq2**2 < 0.000001:
            H2.append(h2)
            M1.append(m1)
            M2.append(m2)
            lambda1, lambda2 = eigenvalues((m1, m2))
            if lambda1>0 and lambda2>0:
                H2_stable.append(h2)
                M1_stable.append(m1)
                M2_stable.append(m2)


fig = plt.figure()
plt.xlabel(r'$h_1=h_2$', fontsize=15)
plt.ylabel(r'$\bar m_1,\ \bar m_2$', fontsize=15)
plt.xlim([h2_l, h2_r])
plt.ylim([-1.02, 1.02])

plt.plot(H2, M1, '*', markersize = 1, label='1st component', alpha=.5, color='darkturquoise')
plt.plot(H2, M2, '*', markersize = 1, label='2nd component', alpha=.5, color='palegreen')

plt.plot(H2_stable, M1_stable, 'o', markersize = 2, label='1st component stable', alpha=.5, color='dodgerblue')
plt.plot(H2_stable, M2_stable, 'o', markersize = 2, label='2nd component stable', alpha=.5, color='limegreen')
plt.legend()

plt.title(title, fontsize=15)
plt.tight_layout()
plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.grid()
fig.set_size_inches(18.2/2, 8.35/1.7)
fig.savefig("../figures/2CW_mags_vs_field.png", dpi=300)
plt.show()
fig.clf()
plt.close()

