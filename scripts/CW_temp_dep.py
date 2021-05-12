########################################################
# This script plots the results of the MC simulations  #
# along with the corresponding theoretical predictions #
# for the temperature dependence of magnetistaion      #
# of the Curie-Weiss model.                            #
########################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
import random

def equation_of_state(m):
    return (h + J*m - 1/(2*beta)*np.log((1+m)/(1-m)))**2

def susceptibility(m):
    return beta*(1-m**2)/(1-beta*J*(1-m**2))

fig, axs = plt.subplots(1, 2, figsize=(18.5/1.4, 10.5/2))

h=0
J=1

T_1 = []
M1 = []
T_step = .001
T_l = T_step
T_r = 1.5

bnds = ((-.99999, .99999),)
for T in np.arange(T_l, T_r, T_step):
    beta = 1/T
    for count in range(0, 10000):
        result = optimize.minimize(equation_of_state, random.uniform(bnds[0][0], bnds[0][1]), bounds = bnds, tol=1e-10)
        if result.fun < .000001:
            T_1.append(T)
            M1.append(result.x[0])

axs[0].plot(T_1, M1, 'o', markersize=2, color='orange', label='Theor')

axs[0].axhline(0, color='black', alpha=.3)
axs[0].axvline(0, color='black', alpha=.3)
axs[0].set_xlabel(r'$T$', fontsize=18)
axs[0].set_ylabel(r'$\bar m$', fontsize=18, labelpad=-5)
axs[0].set_xlim([.2, T_r])
axs[0].set_ylim([-1.01, 1.01])

axs[0].set_title(r'$J=1;\ h=0;\ N=10^5$', fontsize=17)

axs[0].grid()


################################################################
h=.01
J=1

T_ = []
M = []
T_step = .001
T_l = T_step
T_r = 1.5

bnds = ((-.99999, .99999),)
for T in np.arange(T_l, T_r, T_step):
    beta = 1/T
    for count in range(0, 10000):
        result = optimize.minimize(equation_of_state, random.uniform(bnds[0][0], bnds[0][1]), bounds = bnds, tol=1e-10)
        if result.fun < .000001:
            T_.append(T)
            M.append(result.x[0])

axs[1].plot(T_, M, 'o', markersize=2, color='orange', label='Theor')

axs[1].axhline(0, color='black', alpha=.3)
axs[1].axvline(0, color='black', alpha=.3)
axs[1].set_xlim([.2, T_r])
axs[1].set_ylim([-1.01, 1.01])
axs[1].set_xlabel(r'$T$', fontsize=18)
axs[1].set_title(r'$J=1;\ h=0.01;\ N=10^5$', fontsize=17)

axs[1].grid()

############# SIMULATION DATA ###################
file_names = ["../data/CW_model/T_scanning/CW__J_1__H_0__INIT_m_-1__T_START_0.001000__T_FIN_2__T_STEP_0.001000__N_NODES_100000__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv",
              "../data/CW_model/T_scanning/CW__J_1__H_0__INIT_m_1__T_START_0.001000__T_FIN_2__T_STEP_0.001000__N_NODES_100000__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv"
]

M_vs_T = []
T_sim = []
M_avrg = []

for file_name in file_names:
    m_vs_T = np.genfromtxt(file_name, delimiter=',')
    M_vs_T.append(m_vs_T)
    # Finding averages
    m1_avrg = []
    T_sim_ = []
    i = 0
    while i < m_vs_T.shape[0]:
        prev_T = m_vs_T[i][0]
        mag1 = []
        while i < m_vs_T.shape[0] and prev_T == m_vs_T[i][0]:
            mag1.append(m_vs_T[i][1])
            i += 1
        m1_avrg.append(np.array(mag1).mean())
        T_sim_.append(prev_T)
    T_sim.append(T_sim_)
    M_avrg.append(m1_avrg)

axs[0].plot(T_sim[0], M_avrg[0], '*', label='Sim', markersize=1, color='blue')
axs[0].plot(T_sim[1], M_avrg[1], '*', markersize=1, color='blue')
axs[0].legend(loc="lower right", fontsize=16)

file_names = ["../data/CW_model/T_scanning/CW__J_1__H_0.010000__INIT_m_-1__T_START_0.001000__T_FIN_2__T_STEP_0.001000__N_NODES_100000__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv",
              "../data/CW_model/T_scanning/CW__J_1__H_0.010000__INIT_m_1__T_START_0.001000__T_FIN_2__T_STEP_0.001000__N_NODES_100000__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv"
]

M_vs_T = []
T_sim = []
M_avrg = []

for file_name in file_names:
    m_vs_T = np.genfromtxt(file_name, delimiter=',')
    M_vs_T.append(m_vs_T)
    # Finding averages
    m1_avrg = []
    T_sim_ = []
    i = 0
    while i < m_vs_T.shape[0]:
        prev_T = m_vs_T[i][0]
        mag1 = []
        while i < m_vs_T.shape[0] and prev_T == m_vs_T[i][0]:
            mag1.append(m_vs_T[i][1])
            i += 1
        m1_avrg.append(np.array(mag1).mean())
        T_sim_.append(prev_T)
    T_sim.append(T_sim_)
    M_avrg.append(m1_avrg)

axs[1].plot(T_sim[0], M_avrg[0], '*', label='Sim', markersize=1, color='blue')
axs[1].plot(T_sim[1], M_avrg[1], '*', markersize=1, color='blue')
axs[1].legend(loc="lower right", fontsize=16)
plt.tight_layout()
plt.savefig("../figures/CW_temp_dep.png")
plt.show()
