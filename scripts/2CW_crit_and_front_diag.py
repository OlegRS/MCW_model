#######################################################
# This script plots the results of the MC simulations #
# along with the theoretical predictions for          #
# magnetisations of the components as functions of h. #
#######################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
import random

fig, axs = plt.subplots(1,2, figsize=(18.5/1.4, 10.5/2))

gamma1 = .5
gamma2 = .5
J11 = 1
J22 = 0.487664
J12 = -1.615554

h2 = 0.669723

h1 = h1_l = .4
h1_r = .85
h1_step = .00005

def equations_of_state(m):
    eq1 = h1 + J11*gamma1*m[0] + J12*gamma2*m[1] - 1/2*np.log((1+m[0])/(1-m[0]))
    eq2 = h2 + J22*gamma2*m[1] + J12*gamma1*m[0] - 1/2*np.log((1+m[1])/(1-m[1]))
    return (eq1, eq2)

def Jacobian_matrix(m):
    return ((J11*gamma1-1/(1-m[0]**2), J12*gamma2), (J12*gamma1, J22*gamma2-1/(1-m[1]**2)))

# Scanning h1
H1 = []
M1 = []
M2 = []

bnds = ((-.99999, .99999), (-.99999, .99999))
m1_prev, m2_prev = random.uniform(bnds[0][0], bnds[0][1]), random.uniform(bnds[1][0], bnds[1][1])
for count in range(0, 200): # Computing good solver initialisation
    m1, m2 = optimize.fsolve(equations_of_state, [random.uniform(bnds[0][0], bnds[0][1]), random.uniform(bnds[1][0], bnds[1][1])], fprime = Jacobian_matrix)
    eq1, eq2 = equations_of_state([m1, m2])
    if eq1**2 + eq2**2 < 0.000001:
        H1.append(h1)
        M1.append(m1)
        M2.append(m2)
        m1_prev = m1
        m2_prev = m2

for h1 in np.arange(h1_l, h1_r, h1_step):
    completeness = (h1-h1_l)/(h1_r-h1_l)*100
    if completeness % 10 <= h1_step/(h1_r-h1_l)*100:
        print(int(completeness), "% complete")
    for count in range(0, 5):
        m1, m2 = optimize.fsolve(equations_of_state, [m1_prev, m2_prev], fprime = Jacobian_matrix)
        eq1, eq2 = equations_of_state([m1, m2])
        if eq1**2 + eq2**2 < 0.000001:
            H1.append(h1)
            M1.append(m1)
            M2.append(m2)
            m1_prev = m1
            m2_prev = m2

axs[0].set_xlabel(r'$h_1$', fontsize=18)
axs[0].set_ylabel(r'$\bar m_1,\ \bar m_2$', fontsize=18, labelpad=-3)
axs[0].set_xlim([h1_l, h1_r])
axs[0].set_ylim([-1.01, 1.01])
axs[0].plot(H1, M1, label='1st component theor', color='dodgerblue', linewidth=2, alpha=.8)
axs[0].plot(H1, M2, label='2nd component theor', color='limegreen', linewidth=2, alpha=.8)


axs[0].set_title(r"$J_{11}=1,\ J_{22}=0.487664,\ J_{12}=-1.615554;\ h_2=0.669723;\ N_1=N_2=10^5$", fontsize=13)

axs[0].axhline(0, color='black')
axs[0].axvline(0, color='black')

################## SIMULATION DATA #######################
file_names = ["../data/2CW_model/h1_scanning/2CW__J111__J22_0.487664__J12_-1.615554__H1_START_0__H1_FIN_1.500000__H1_STEP_0.001000__H2_0.669723__N1_100000__N2_100000__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv"]

M_vs_H = []
H = []
M1_avrg = []
M2_avrg = []

for file_name in file_names:
    m_vs_h = np.genfromtxt(file_name, delimiter=',')
    M_vs_H.append(m_vs_h)
    # Finding averages
    m1_avrg = []
    m2_avrg = []
    h = []
    i = 0
    while i < m_vs_h.shape[0]:
        prev_h = m_vs_h[i][0]
        mag1 = []
        mag2 = []
        while i < m_vs_h.shape[0] and prev_h == m_vs_h[i][0]:
            mag1.append(m_vs_h[i][1])
            mag2.append(m_vs_h[i][2])
            i += 1
        m1_avrg.append(np.array(mag1).mean())
        m2_avrg.append(np.array(mag2).mean())
        h.append(prev_h)

axs[0].plot(h, m1_avrg, 'o', label='1st component sim', markersize=3, color='darkblue', alpha=0.4, zorder=15)
axs[0].plot(h, m2_avrg, '*', label='2nd component sim', markersize=3.7, color='darkgreen', alpha=0.4, zorder=15)
axs[0].plot(0.563408, .3, '+', markersize=14, color='violet', markeredgewidth=1.5, zorder=10, alpha=.9,label='Critical point')
axs[0].plot(0.563408, .5, '+', markersize=14, color='violet', markeredgewidth=1.5, zorder=10, alpha=.9)


axs[0].legend(loc="lower right", fontsize=14)
axs[0].grid()


#############################################
################ TWO WAYS SCAN ##############
#############################################
gamma1 = .5
gamma2 = 1-gamma1

J11 = 1
J22 = .487664
J12 = -1.615554

h2_l = -1
h2_r = 1
h2_step = .001


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
    h1 = h2
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


axs[1].set_xlabel(r'$h_1=h_2$', fontsize=18)
axs[1].set_xlim([h2_l, h2_r])
axs[1].set_ylim([-1.01, 1.01])

axs[1].plot(H2, M1, '*', markersize = 1, label='1st component unstable', alpha=.5, color='darkturquoise')
axs[1].plot(H2, M2, '*', markersize = 1, label='2nd component unstable', alpha=.5, color='palegreen')

axs[1].plot(H2_stable, M1_stable, 'o', markersize = 2, label='1st component stable', alpha=.5, color='dodgerblue')
axs[1].plot(H2_stable, M2_stable, 'o', markersize = 2, label='2nd component stable', alpha=.5, color='limegreen')

################## SIMULATION DATA #######################
file_names = ["../data/2CW_model/h1=h2_scanning/2CW__J111__J22_0.487664__J12_-1.615554__H1_H2_START_-1.100000__H1_H2_FIN_1.100000__H1_H2_STEP_0.010000__N1_100000__N2_100000__ITERS_PER_NODE_10__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv",
              "../data/2CW_model/h1=h2_scanning/2CW__J111__J22_0.487664__J12_-1.615554__H1_H2_START_1.100000__H1_H2_FIN_-1.100000__H1_H2_STEP_-0.010000__N1_100000__N2_100000__ITERS_PER_NODE_10__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv"
]

M_vs_H = []
H = []
M1_avrg = []
M2_avrg = []

for file_name in file_names:
    print(file_name)
    m_vs_h = np.genfromtxt(file_name, delimiter=',')
    M_vs_H.append(m_vs_h)
    # Finding averages
    m1_avrg = []
    m2_avrg = []
    h = []
    i = 0
    while i < m_vs_h.shape[0]:
        prev_h = m_vs_h[i][0]
        mag1 = []
        mag2 = []
        while i < m_vs_h.shape[0] and prev_h == m_vs_h[i][0]:
            mag1.append(m_vs_h[i][1])
            mag2.append(m_vs_h[i][2])
            i += 1
        m1_avrg.append(np.array(mag1).mean())
        m2_avrg.append(np.array(mag2).mean())
        h.append(prev_h)
    axs[1].plot(h, m1_avrg, 'o', markersize=3, color='darkblue', alpha=0.5)
    axs[1].plot(h, m2_avrg, '*', markersize=3.7, color='darkgreen', alpha=0.5)

axs[1].set_title(r"$J_{11}=1,\ J_{22}=0.487664,\ J_{12}=-1.615554;\ N_1=N_2=10^5$", fontsize=13)

axs[1].axhline(0, color='black')
axs[1].axvline(0, color='black')
axs[1].grid()

plt.tight_layout()
fig.savefig("../figures/2CW_crit_and_front_diag.png", dpi=300)
plt.show()

fig.clf()
plt.close()
