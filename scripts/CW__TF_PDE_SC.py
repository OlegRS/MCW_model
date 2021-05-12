###########################################################
# This script plots the results of the ensemble averaging #
# MC simulations for the Curie-Weiss model together with  #
# the corresponding theoretical predictions. The exact    #
# evaluation of expectation is done by solving the        #
# Burgers' equation using central differences in space    #
# with RK4 time stepping.                                 #
###########################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

N = 50
h_min = -.1
h_max = .1

fig, ax = plt.subplots(figsize=(5*2, 3*2))

##################### RK4 PDE SOLVER ############################
x_min = -5
x_max = 5
nx = 10000
dx = (x_max-x_min)/nx
t2 = 2
nJ = 50000
dJ = t2/nJ
t3 = 1
nt3 = 100000
dt3 = t3/nt3
x = np.arange(x_min, x_max, dx)

nu = 1/(2*N)

# Setting the initial conditions
u = np.tanh(x)

### Defining auxiliary functions ###
def d1x(y):
    D1x = np.zeros(nx)
    D1x[1:nx-1]=(y[2:nx]-y[0:nx-2])/(2*dx)
    D1x[0]=0
    D1x[nx-1]=0
    return D1x

def F1(u_): # Right hand side of the 1st equation of the hierarchy
    return d1x(.5*u_**2 + nu*d1x(u_))

def RK4(u_, F, dt):
    K1 = F(u_)
    K2 = F(u_+.5*dt*K1)
    K3 = F(u_+.5*dt*K2)
    K4 = F(u_+dt*K3)
    return u+dt/6.*(K1 + 2*K2 + 2*K3 + K4)

### RK4 time stepping for 1st eqn of hierarchy ###
for j in range(0, nJ):
    u = RK4(u, F1, dJ)
##################### RK4 PDE SOLVER END ##########################
file_names = ["../data/CW_model/h_scanning/CW__J_2__H_START_-0.100000__H_FIN_-0.080000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_-0.080000__H_FIN_-0.060000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_-0.060000__H_FIN_-0.040000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_-0.040000__H_FIN_-0.020000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_-0.020000__H_FIN_0__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_0__H_FIN_0.020000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_0.020000__H_FIN_0.040000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_0.040000__H_FIN_0.060000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_0.060000__H_FIN_0.080000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv",
              "../data/CW_model/h_scanning/CW__J_2__H_START_0.080000__H_FIN_0.100000__H_STEP_0.001000__N_NODES_50__ITERS_PER_NODE_1000__N_AVRGING_10000__RANDOM_INITIALISATION_0__N_DYNAMIC_NODES_50.csv"
]

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

    H.extend(h)
    M1_avrg.extend(m1_avrg)

plt.plot(H, M1_avrg, 'o', label='Simulation', markersize=3, color='blue', zorder=10, alpha=.6)

####################### SC CURVE #####################################
T = 1.
J=2
# Plotting self-consistency equation
delta = 0.0001
Frange = np.arange(h_min, h_max+delta, delta)
mrange = np.arange(-1.1, 1.1, delta)
M, F_ = np.meshgrid(mrange, Frange)
G = M - np.tanh(J/T*M + F_/T)
CS = plt.contour(F_, M, G, [0], colors='cyan')
CS.collections[0].set_label('Self-consistency equation')
####################### SC CURVE END #################################

#################### TRANSITION FORMULA ####################
h = np.arange(h_min, h_max+delta, delta)
m_TF = np.tanh(N*h)
plt.plot(h, m_TF, color='green', label='Transition formula')
################### TRANSITION FORMULA END ##################
plt.plot(x, u, color='red', label='Exact solution', alpha=1)

plt.axvline(0, color='black')
plt.axhline(0, color='black')
plt.grid()
lgnd = plt.legend(loc=[.02, .67], prop={'size': 16})
plt.title(r"$N=50,\ J=2$", size=20)

plt.xlim([h_min, h_max])
plt.ylim([-1.01, 1.01])
plt.xlabel(r'$h$', fontsize=20)
plt.ylabel(r'$\bar m$', fontsize=20)

##### ZOOM ######
axins = inset_axes(ax, width="50%", height="50%", bbox_to_anchor=(.01, .05, .5, .82), bbox_transform=ax.transAxes, loc="lower left")

CS = axins.contour(F_, M, G, [0], colors='cyan')
axins.plot(x, u, color='red', alpha=.9, zorder=10)
axins.plot(h, m_TF, color='green')
axins.plot(H, M1_avrg, 'o', markersize=3, color='blue')

x1, x2, y1, y2 = -.014, -.003, -.55, -.15 # specify the limits
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.yticks(visible=False)
plt.xticks(visible=False)

mark_inset(ax, axins, loc1=1, loc2=4, fc="none", ec="0.5")

axins = inset_axes(ax, width="50%", height="50%", bbox_to_anchor=(.64, .51, .7, .75), bbox_transform=ax.transAxes, loc="lower left")

CS = axins.contour(F_, M, G, [0], colors='cyan')
axins.plot(x, u, color='red', alpha=.9, zorder=10)
axins.plot(h, m_TF, color='green')
axins.plot(H, M1_avrg, 'o', markersize=3, color='blue')

x1, x2, y1, y2 = .052, .09, .945, 1 # specify the limits
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.yticks(visible=False)
plt.xticks(visible=False)

mark_inset(ax, axins, loc1=1, loc2=2, fc="none", ec="0.5")

plt.tight_layout()
plt.savefig('../figures/CW__TF_PDE_SC.png', dpi=300)
plt.show()
fig.clear()
plt.close(fig)
