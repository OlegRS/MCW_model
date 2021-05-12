###################################################
# This script plots the results of time averaging #
# MC simulations for the Curie-Weiss model along  #
# with the corresponding theoretical predictions. #
###################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np

file_names = ["../data/CW_model/h_scanning/CW__J_2__H_START_-1.530000__H_FIN_1.530000__H_STEP_0.010000__N_NODES_10000__ITERS_PER_NODE_10__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_10000.csv"]

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 3.2))

M_vs_H = []
H = []
M1_avrg = []

for file_name in file_names:
    m_vs_h = np.genfromtxt(file_name, delimiter=',')
    M_vs_H.append(m_vs_h)
    # Finding averages
    m1_avrg = []
    h = []
    i = 0
    while i < m_vs_h.shape[0]:
        prev_h = m_vs_h[i][0]
        mag1 = []
        while i < m_vs_h.shape[0] and prev_h == m_vs_h[i][0]:
            mag1.append(m_vs_h[i][1])
            i += 1
        m1_avrg.append(np.array(mag1).mean())
        h.append(prev_h)
    H.append(h)
    M1_avrg.append(m1_avrg)
    
# Plotting
axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')
axs[0].set_xlabel(r'$h$', fontsize=15)
axs[0].set_ylabel(r'$\bar m$', fontsize=15, labelpad=-5)
axs[0].set_xlim([-1.5, 1.5])
axs[0].set_ylim([-1.03, 1.03])

axs[0].plot(H[0], M1_avrg[0], 'o', label='Simulation', markersize=3, color='blue', zorder=100, alpha=.5)

####################### THEORETICAL CURVE #####################################
T = 1.
# Plotting self-consistency equation
axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')

J=2
delta = 0.0025
Frange = np.arange(-1.5, 1.5, delta)
mrange = np.arange(-1.1, 1.1, delta)
M, F = np.meshgrid(mrange, Frange)
G = M - np.tanh(J/T*M + F/T)
CS = axs[0].contour(F, M, G, [0], colors='cyan', zorder=50)
CS.collections[0].set_label('SC equation')
CS.collections[0].set_linewidth(2)
####################### THEORETICAL CURVE END #################################
axs[0].grid()

###################### BACK SCAN ########################
file_names = ["../data/CW_model/h_scanning/CW__J_2__H_START_1.530000__H_FIN_-1.530000__H_STEP_-0.010000__N_NODES_10000__ITERS_PER_NODE_10__N_AVRGING_10__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_10000.csv"]

M_vs_H = []
H = []
M1_avrg = []

for file_name in file_names:
    m_vs_h = np.genfromtxt(file_name, delimiter=',')
    M_vs_H.append(m_vs_h)

    # Finding averages
    m1_avrg = []
    
    h = []
    i = 0
    while i < m_vs_h.shape[0]:
        prev_h = m_vs_h[i][0]
        mag1 = []
        while i < m_vs_h.shape[0] and prev_h == m_vs_h[i][0]:
            mag1.append(m_vs_h[i][1])
            i += 1
        m1_avrg.append(np.array(mag1).mean())

        h.append(prev_h)

    H.append(h)
    M1_avrg.append(m1_avrg)
    
# Plotting
# fig = plt.figure(dpi=200)

axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')

axs[1].set_xlabel('$h$', fontsize=15)
# axs[1].set_ylabel(r'$\langle m\rangle$', fontsize=15)

axs[1].set_xlim([-1.5, 1.5])
axs[1].set_ylim([-1.03, 1.03])


####################### THEORETICAL CURVE #####################################

T = 1.
# Plotting self-consistency equation
axs[1].axvline(0, color='black')
axs[1].axhline(0, color='black')

CS = axs[1].contour(F, M, G, [0], colors='cyan', zorder=50)
CS.collections[0].set_label('SC equation')
CS.collections[0].set_linewidth(2)

axs[1].plot(H[0], M1_avrg[0], 'o', label='Simulation', markersize=3, color='blue', zorder=100, alpha=.5)

##### Exact thermodynamic solution #####
x_min_, x_max_ = -.1, .1
H = np.arange(x_min_, x_max_, .0001)
m_MIN = []
for h in H:
    def SC_equation(m):
        return m - np.tanh(J/T*m+ h/T)
    m_ = np.array([fsolve(SC_equation, -1)[0], fsolve(SC_equation, 1)[0]])
    def F(m):
        return -J*m**2/2 - h*m + .5*T*np.log((1+m)**(1+m)*(1-m)**(1-m))
    m0 = m_[~(np.triu(np.abs(m_[:,None] - m_) < .001, 1)).any(0)] # Distinct solution of SC equation
    m_ = np.array([l for l in m0 if (l<1 and l>-1)])
    m_min = m_[list(F(m_)).index(min(F(m_)))]
    m_MIN.append(m_min)

axs[0].plot(H, m_MIN, color="red", label='Phase transition', zorder=10, linewidth=2)
axs[1].plot(H, m_MIN, color="red", label='Phase transition', zorder=10, linewidth=2)

####################### THEORETICAL CURVE END #################################
axs[1].grid()
axs[0].legend(loc=[.54,.62], prop={'size': 11})
axs[1].legend(loc=[.54,.62], prop={'size': 11})
# axs[0].set_title(r'$J=2;\ T=1;\ N=10^4$')
# axs[1].set_title(r'$J=2;\ T=1;\ N=10^4$')

plt.tight_layout()
plt.savefig('../figures/CW_2_way_scan.png', dpi=300)
plt.show()
