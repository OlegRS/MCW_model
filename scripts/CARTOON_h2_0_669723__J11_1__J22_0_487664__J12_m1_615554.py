###############################################################
# This script produces frames from the animation available at #
# https://youtu.be/2gJrtE2QUaY                                #
###############################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import time
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
import random

gamma1 = .5
gamma2 = 1-gamma1

J11 = 1
J22 = .487664
J12 = -1.615554

delta = 0.001
m1 = np.arange(-.99999, .99999, delta)
m2 = np.arange(-.99999, .99999, delta)

h1_min, h1_max = -1, 1
h2_min, h2_max = -1, 1
h2_fixed = 0.669723

NUM_FRAMES = 300
h2_scan_min = h1_min
h2_scan_max = h1_min + (h1_max - h1_min)/NUM_FRAMES

H_curve_plot = []
M_curve_plot = []
for frame in range(0, NUM_FRAMES):
    print("Computing frame ", frame, "out of ", NUM_FRAMES)
################## FINDING NUMERICAL SOLUTIONS ################################
    # Defining h-curve
    H_curve = []
    for h1 in np.arange(h2_scan_min, h2_scan_max, 0.001):
        H_curve.append([h1, h2_fixed])
        H_curve_plot.append([h1, h2_fixed])
    H_curve = np.array(H_curve)
    H_curve_plot_array = np.array(H_curve_plot)
    
    def equation_of_state(m):
        eq1 = h1 + J11*gamma1*m[0] + J12*gamma2*m[1] - 1/2*np.log((1+m[0])/(1-m[0]))
        eq2 = h2 + J22*gamma2*m[1] + J12*gamma1*m[0] - 1/2*np.log((1+m[1])/(1-m[1]))
        return eq1**2 + eq2**2

    M_curve = []
    bnds = ((-.99999, .99999), (-.99999, .99999))
    for h1, h2 in H_curve:
        for count in range(0, 200): #200
            result = optimize.minimize(equation_of_state, [random.uniform(bnds[0][0], bnds[0][1]), random.uniform(bnds[1][0], bnds[1][1])], bounds = bnds, tol=1e-10)
            if result.fun < .000001:
                M_curve.append([result.x[0], result.x[1]])
                M_curve_plot.append([result.x[0], result.x[1]])

################### FINDING NUMERICAL SOLUTIONS END ############################

    M1, M2 = np.meshgrid(m1, m2)

    J11_ = J11-1/(gamma1*(1-M1**2))
    J22_ = J22-1/(gamma2*(1-M2**2))
    lambda1 = -(J11_ + J22_ + np.sqrt((J11_-J22_)**2 + 4*J12**2))/2
    lambda2 = -(J11_ + J22_ - np.sqrt((J11_-J22_)**2 + 4*J12**2))/2
    DET = lambda1*lambda2

    lambda1_positive = lambda1 > 0
    lambda2_positive = lambda2 > 0

    F_SHAPE = np.empty((M1.shape[0], M1.shape[1]))

    for i in range(M1.shape[0]):
        for j in range(M1.shape[1]):
            if lambda1[i, j] > 0 and lambda2[i, j] > 0:
                F_SHAPE[i, j] = 1
            if (lambda1[i, j] > 0 and lambda2[i, j] < 0) or (lambda1[i, j] < 0 and lambda2[i, j] > 0):
                F_SHAPE[i, j] = 0
            if lambda1[i, j] < 0 and lambda2[i, j] < 0:
                F_SHAPE[i, j] = -1

    H1 = -gamma1*J11*M1 - gamma2*M2*J12 + np.arctanh(M1)
    H2 = -gamma2*J22*M2 - gamma1*M1*J12 + np.arctanh(M2)

    fig ,axs = plt.subplots(nrows=1, ncols=2, figsize=(18.2/1.7, 8.35/1.7))
    
    axs[1].axis('equal')
    axs[1].axhline(0, linewidth=.5, color='black')
    axs[1].axvline(0, linewidth=.5, color='black')
    cs = axs[1].contourf(M1, M2, F_SHAPE)
    axs[1].plot([M_curve_plot[i][0] for i in range(0, M_curve_plot.__len__())], [M_curve_plot[i][1] for i in range(0, M_curve_plot.__len__())], '*', markersize=2, color='cyan')
    axs[1].plot([M_curve[i][0] for i in range(0, M_curve.__len__())], [M_curve[i][1] for i in range(0, M_curve.__len__())], '*', markersize=3, color='blue', zorder=110)
    axs[1].contour(M1, M2, DET, 0, colors='r')
    axs[1].set_xlim([-1, 1])
    axs[1].set_ylim([-1, 1])
    axs[1].set_xlabel(r'$\bar m_1$', fontsize=15)
    axs[1].set_ylabel(r'$\bar m_2$', fontsize=15, labelpad=-8)


    axs[1].text(-.32, 0, 'Saddle point\n(UNSTABLE)', color='red', zorder=100, fontsize=12, weight='bold', alpha=.8)
    axs[1].text(-.8, -.85, 'Local minimum\n    (STABLE)', color='red', zorder=100, fontsize=12, weight='bold', alpha=.8)

    axs[1].grid()

    axs[0].axis('equal')
    axs[0].axhline(0, linewidth=1, color='black')
    axs[0].axvline(0, linewidth=1, color='black')
    axs[0].plot(H_curve_plot_array[:, 0], H_curve_plot_array[:, 1], linewidth=2, color='cyan')
    axs[0].plot(H_curve[:, 0], H_curve[:, 1], linewidth=3, color='blue', zorder=110)
    axs[0].contour(H1, H2, DET, 0, colors='r')
    axs[0].set_xlim([-1, 1])
    axs[0].set_ylim([h2_min, h2_max])
    axs[0].set_xlabel('$h_1$', fontsize=15)
    axs[0].set_ylabel('$h_2$', fontsize=15, labelpad=-8)

    axs[0].text(.01, .05, '3-2', color='blue', zorder=100, fontsize=11)
    axs[0].text(-.5, .3, '1-1', color='blue', zorder=100, fontsize=11)

    axs[0].grid()

    fig.suptitle("$J_{11}=$"+str(J11)+"$,\ J_{22}=$"+str(J22)+"$,\ J_{12}=$"+str(J12)+"$;\ \ \gamma_1=$"+str(gamma1)+"$,\ \gamma_2=$"+str(gamma2) + "$;\ h_1=$" + '{:<06}'.format(round(h2_scan_max,3)) + "$,\ h_2=$" + str(h2_fixed), fontsize=13, y=0.94)
    # plt.show()

    fig.savefig("../figures/CARTOON_FRAMES/crit_h2_const/"+time.strftime("%Y%m%d-%H%M%S")+"J11_"+str(J11)+"_J22_"+str(J22)+"_J12_"+str(J12)+"_gamma1_"+str(gamma1)+"_gamma2_"+str(gamma2)+"_h1_min_"+str(h1_min)+"_h1_max_"+str(h1_max)+"_h2_min_"+str(h2_min)+"_h2_max_"+str(h2_max)+"_delta_"+str(delta)+".png", dpi=250)

    fig.clf()
    plt.close()

    h2_scan_min = h2_scan_max
    h2_scan_max += (h1_max - h1_min)/NUM_FRAMES
