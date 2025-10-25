#           Program "plane_wave_trans.py"
#
# Comparison of the analytical expression and the expansion of the primary
# incident plane wave after translation from the laboratory frame into the 
# body frame (rotation - shift - rotation). Calculation is performed along 
# the surface of a shifted sphere with radius "a1", and for "theta_1" in  
# [0, pi] (in steps of 10 degree)!
#

print()
print()
print("       Translation of the primary incident plane wave!")
print()
print()

import numpy as np
import basics as bas
import scipy as scp
import scipy.special as scs


# Input parameters:
    
a1 = float(input(' radius of the shifted sphere ... a1 [mm]: '))
b = float(input(' shift along the new z-direction after 1. rotation ... b [mm]: '))
alpha = float(input(' Eulerian angle ... alpha [deg.]: '))
theta_p = float(input(' Eulerian angle ... theta_p [deg.]: '))
ncut = int(input(" truncation parameter ...  ncut: "))
print()
print()

k0 = 1.0
k0b = k0 * b
beta_1 = k0 * a1
zi = 0. + 1.0j
pf = - zi / k0

theta_l = 180.0
theta_l1 = 19
theta = np.linspace(0.0, theta_l, theta_l1)
r_theta = theta * scp.pi / 180.
r_thetap = theta_p * scp.pi / 180.


# Calculation of the separation matrix:

sm_l_i_j = []   
for l in range(0, ncut + 1):   
    sm_i_j = []
    for i in range(0, ncut + 1):
        sm_j = []
        for j in range(0, ncut + 1):
            sm_j = sm_j + [bas.SM(l, i, j, k0b)]
        sm_i_j = sm_i_j + [sm_j]
    sm_l_i_j = sm_l_i_j + [sm_i_j]
sm_l_i_j_real = np.real(sm_l_i_j)


# Calculation of the primary incident plane wave in the body frame:

for i in range(0, theta_l1):
    sum_zw = 0.0
    ana_zw = np.exp(zi * k0b * np.cos(r_thetap) + zi * beta_1 *\
                    np.cos(r_theta[i]))
    for n_zs in range(0, ncut + 1):
        bes = scs.spherical_jn(n_zs,beta_1)
        for l1 in range(-n_zs, n_zs + 1):
            Y_sh = scs.sph_harm(l1,n_zs,0.,r_theta[i])
            for ls in range(-ncut, ncut + 1):
                s_ls = sm_l_i_j_real[np.abs(ls)]
                s_ls_nzs = s_ls[n_zs]
                for n in range(np.abs(ls), ncut + 1):
                    d_n = np.sqrt(4. * np.pi * (2 * n + 1)) * zi**n
                    d_n = (-1)**n * d_n
                    Dm_1 = bas.drm(n, ls, 0, 0., -theta_p, -alpha)
                    Dm_2 = bas.drm(n_zs, l1, ls, alpha, theta_p, 0.)
                    S_m = s_ls_nzs[n]
                    sum_zw = sum_zw + d_n * Dm_1 * Dm_2 * S_m * bes * Y_sh                 
    print(' real ana vs. sum:')    
    print(' theta = ',theta[i],'deg.: ',np.real(ana_zw),np.real(sum_zw))
    print(' imag ana vs. sum')
    print(' theta = ',theta[i],'deg.: ',np.imag(ana_zw),np.imag(sum_zw))
    print()
print(' Press Enter to finish')
input()
