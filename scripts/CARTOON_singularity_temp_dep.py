######################################################################
# This script produces frames from the animation available at        #
# https://youtu.be/amoPQkPNDdI. Colours on the first 91 frames       #
# are corrected using                                                #
# ../figures/CARTOON_FRAMES/singularity_temp_dep/color_correction.sh #
# to make (meta)stability region appear consistently yellow. The     #
# colour correction script requires ImageMagic.                      #
######################################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import time
import matplotlib.pyplot as plt
import numpy as np

gamma1 = .5
gamma2 = .5
J_11 = 3
J_22 = 3
J_12 = -.42

T_min = 0
T_max = 1.5

delta = .0005
m1 = np.arange(-.99999, 1, delta)
m2 = np.arange(-.99999, 1, delta)

h1_min = -1
h1_max = 1
h2_min = -1
h2_max = 1

M1, M2 = np.meshgrid(m1, m2)

for T in np.arange(1.8, 0, -.001):
    J11 = J_11/T
    J22 = J_22/T
    J12 = J_12/T
    J11_ = J11-1/(gamma1*(1-M1**2))
    J22_ = J22-1/(gamma2*(1-M2**2))
    lambda1 = -(J11_ + J22_ + np.sqrt((J11_-J22_)**2 + 4*J12**2))/2
    lambda2 = -(J11_ + J22_ - np.sqrt((J11_-J22_)**2 + 4*J12**2))/2
    DET = lambda1*lambda2

    # DET_positive = DET > 0
    lambda1_positive = lambda1 > 0
    lambda2_positive = lambda2 > 0


    F_SHAPE = np.empty((M1.shape[0], M1.shape[1]))

    for i in range(M1.shape[0]):
        for j in range(M1.shape[1]):
            if lambda1[i, j] > 0 and lambda2[i, j] > 0:
                F_SHAPE[i, j] = 1
            if (lambda1[i, j] > 0 and lambda2[i, j] < 0) or (lambda1[i, j] < 0 and lambda2[i, j] > 0):
                F_SHAPE[i, j] = -1
            if lambda1[i, j] < 0 and lambda2[i, j] < 0:
                F_SHAPE[i, j] = 0

    H1 = -gamma1*J11*M1 - gamma2*M2*J12 + np.arctanh(M1)
    H2 = -gamma2*J22*M2 - gamma1*M1*J12 + np.arctanh(M2)

    fig ,axs = plt.subplots(nrows=1, ncols=2, figsize=(18.2/1.7, 8.35/1.7))
    axs[1].axis('equal')
    axs[1].axhline(0, linewidth=.5, color='black')
    axs[1].axvline(0, linewidth=.5, color='black')

    cs = axs[1].contourf(M1, M2, F_SHAPE)

    axs[1].contour(M1, M2, DET, 0, colors='r')
    axs[1].set_xlim([-1, 1])
    axs[1].set_ylim([-1, 1])
    axs[1].set_xlabel(r'$\bar m_1$', fontsize=15)
    axs[1].set_ylabel(r'$\bar m_2$', fontsize=15, labelpad=-8)
    axs[1].grid()

    axs[0].axis('equal')
    axs[0].axhline(0, linewidth=1, color='black')
    axs[0].axvline(0, linewidth=1, color='black')
    axs[0].contour(H1, H2, DET, 0, colors='r')
    axs[0].set_xlim([h1_min, h1_max])
    axs[0].set_ylim([h2_min, h2_max])
    
    axs[0].set_xlabel(r'$h_1$', fontsize=15)
    axs[0].set_ylabel(r'$h_2$', fontsize=15, labelpad=-8)
    axs[0].grid()

    fig.suptitle(r'$J_{11}=$'+str(J_11)+r'$,\ J_{22}=$'+str(J_22)+r'$,\ J_{12}=$'+str(J_12)+r'$;\ T=$'+'{:<05}'.format(round(T,4))+r'$;\ \gamma_1=$'+str(gamma1)+r'$,\ \gamma_2=$'+str(gamma2), fontsize=13, y=0.94)

    # plt.show()

    fig.savefig("../figures/CARTOON_FRAMES/singularity_temp_dep/"+time.strftime("%Y%m%d-%H%M%S")+"J11_"+str(J11)+"_J22_"+str(J22)+"_J12_"+str(J12)+"_gamma1_"+str(gamma1)+"_gamma2_"+str(gamma2)+"_h1_min_"+str(h1_min)+"_h1_max_"+str(h1_max)+"_h2_min_"+str(h2_min)+"_h2_max_"+str(h2_max) + "_T_min_"+str(T_min)+"_T_max_"+str(T_max) +"_delta_"+str(delta)+".png", dpi=300)

    fig.clf()
    plt.close()
