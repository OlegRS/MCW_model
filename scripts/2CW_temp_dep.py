#########################################################
# This script plots temperature dependence of the full  #
# specific magnetisation of the 2-component Curie-Weiss #
# model together with that for gamma1*m1-gamma2*m2.     #
#########################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
import random

gamma1 = .5
gamma2 = 1 - gamma1
T = 1
T_r = 2.2
T_l = 0.001
T_step = -0.001 
J11 = 3
J22 = 3
J12 = -.42

h = 0
h1 = h
h2 = h

fig, axs = plt.subplots(1, 2, figsize=(18.5/1.4, 10.5/2))

def equation_of_state(m):
    eq1 = (h1 + J11*gamma1*m[0] + J12*gamma2*m[1])/T - 1/2*np.log((1+m[0])/(1-m[0]))
    eq2 = (h2 + J22*gamma2*m[1] + J12*gamma1*m[0])/T - 1/2*np.log((1+m[1])/(1-m[1]))
    return eq1**2 + eq2**2

def eigenvalues(m):
    J11_ = J11/T-1/(gamma1*(1-m[0]**2))
    J22_ = J22/T-1/(gamma2*(1-m[1]**2))
    J12_ = J12/T
    lambda1 = -(J11_ + J22_ + np.sqrt((J11_-J22_)**2 + 4*J12_**2))/2
    lambda2 = -(J11_ + J22_ - np.sqrt((J11_-J22_)**2 + 4*J12_**2))/2
    return lambda1, lambda2


T_ = []
M1 = []
M2 = []
M_full = []
M_diff = []
T_stable = []
M1_stable = []
M2_stable = []
M_full_stable = []
M_diff_stable = []

bnds = ((-.9999999999, .9999999999), (-.9999999999, .9999999999))
# Scanning J12
for T in np.arange(T_r, T_l, T_step):
    for count in range(0, 10000):
        result = optimize.minimize(equation_of_state, [random.uniform(bnds[0][0], bnds[0][1]), random.uniform(bnds[1][0], bnds[1][1])], bounds = bnds, tol=1e-10)
        if result.fun < .000001:
            T_.append(T)
            M1.append(result.x[0])
            M2.append(result.x[1])
            M_full.append(gamma1*result.x[0] + gamma2*result.x[1])
            M_diff.append(gamma1*result.x[0] - gamma2*result.x[1])
            lambda1, lambda2 = eigenvalues((result.x[0], result.x[1]))
            if lambda1>0 and lambda2>0:
                T_stable.append(T)
                M1_stable.append(result.x[0])
                M2_stable.append(result.x[1])
                M_full_stable.append(gamma1*result.x[0] + gamma2*result.x[1])
                M_diff_stable.append(gamma1*result.x[0] - gamma2*result.x[1])

################## SIMULATION DATA #######################

file_names = ["../data/2CW_model/T_scanning/2CW__J113__J22_3__J12_-0.420000__T_START_0.001000__T_FIN_2.500000__T_STEP_0.001000__H1_0__H2_0__N1_100000__N2_100000__INIT_m1_-1__INIT_m2_-1__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv",
              "../data/2CW_model/T_scanning/2CW__J113__J22_3__J12_-0.420000__T_START_0.001000__T_FIN_2.500000__T_STEP_0.001000__H1_0__H2_0__N1_100000__N2_100000__INIT_m1_-1__INIT_m2_1__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv",
              "../data/2CW_model/T_scanning/2CW__J113__J22_3__J12_-0.420000__T_START_0.001000__T_FIN_2.500000__T_STEP_0.001000__H1_0__H2_0__N1_100000__N2_100000__INIT_m1_1__INIT_m2_-1__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv",
              "../data/2CW_model/T_scanning/2CW__J113__J22_3__J12_-0.420000__T_START_0.001000__T_FIN_2.500000__T_STEP_0.001000__H1_0__H2_0__N1_100000__N2_100000__INIT_m1_1__INIT_m2_1__ITERS_PER_NODE_500__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_1.csv"
]

M_vs_T = []

legend_set = False
for file_name in file_names:
    m_vs_T = np.genfromtxt(file_name, delimiter=',')
    M_vs_T.append(m_vs_T)
    # Finding averages
    m1_avrg = []
    m2_avrg = []
    m_tot_avrg = []
    m_diff_avrg = []
    T = []
    i = 0
    while i < m_vs_T.shape[0]:
        prev_T = m_vs_T[i][0]
        mag1 = []
        mag2 = []
        m_tot = []
        m_diff = []
        while i < m_vs_T.shape[0] and prev_T == m_vs_T[i][0]:
            mag1.append(m_vs_T[i][1])
            mag2.append(m_vs_T[i][2])
            m_tot.append(gamma1*m_vs_T[i][1] + gamma2*m_vs_T[i][2])
            m_diff.append(gamma1*m_vs_T[i][1] - gamma2*m_vs_T[i][2])
            i += 1
        m1_avrg.append(np.array(mag1).mean())
        m2_avrg.append(np.array(mag2).mean())
        m_tot_avrg.append(np.array(m_tot).mean())
        m_diff_avrg.append(np.array(m_diff).mean())
        T.append(prev_T)
    if not legend_set:
        axs[0].plot(T, m_tot_avrg, '*', markersize=1, color='blue', alpha=1, zorder=100, label='Sim')
        axs[1].plot(T, m_diff_avrg, '*', markersize=1, color='blue', alpha=1, zorder=100, label='Sim')
        legend_set = True
    else:
        axs[0].plot(T, m_tot_avrg, '*', markersize=1, color='blue', alpha=1, zorder=100)
        axs[1].plot(T, m_diff_avrg, '*', markersize=1, color='blue', alpha=1, zorder=100)

##############################################################################################################

axs[0].axhline(0, color='black')
axs[0].axvline(0, color='black')
axs[0].set_xlabel('$T$', fontsize=18)
axs[0].set_ylabel(r'$\gamma_1\bar m_1 + \gamma_2\bar m_2$', fontsize=18)
axs[0].set_xlim([.3, T_r])
axs[0].set_ylim([-1.01, 1.01])
axs[0].plot(T_, M_full, 'o', markersize=2, label='Theor', color='orange')
axs[0].legend(loc="upper right", fontsize=16)
axs[0].grid()


axs[1].axhline(0, color='black')
axs[1].axvline(0, color='black')
axs[1].set_xlabel('$T$', fontsize=18)
axs[1].set_ylabel(r'$\gamma_1\bar m_1 - \gamma_2\bar m_2$', fontsize=18)
axs[1].set_xlim([.3, T_r])
axs[1].set_ylim([-1.01, 1.01])
axs[1].plot(T_, M_diff, 'o', markersize=2, label='Theor', color='orange')
axs[1].legend(loc="upper right", fontsize=16)
axs[1].grid()

# plt.suptitle(r'$J_{11}=$'+str(J11)+r'$,\ J_{22}=$'+str(J22)+r'$,\ J_{12}=$'+str(J12)+r';$\ h_1=$'+str(h1)+r'$,\ h_2=$'+str(h2) + r'$;\ N_1=10^5$' + r'$;\ N_2=10^5$', fontsize=17, y=1)

plt.tight_layout()
plt.subplots_adjust(top=.946)

fig.savefig("../figures/2CW_temp_dep.png", dpi=300)

plt.show()
fig.clf()
plt.close()
