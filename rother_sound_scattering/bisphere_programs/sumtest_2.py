#            Program "sumtest_2.py"
#
# This program provides a test of Eq. (2.74) of the manuscript.
#

print()
print()
print("         Test of Eq. (2.74) of the manuscript!")
print()
print()

import numpy as np
import scipy as scp
import scipy.special as scs
import matplotlib.pyplot as plt
import basics as bas
    

# Input parameter:

para = str(input('soft (s) or hard (h) or penetrable (p)? '))

if para == 'p':
    k_p_k = float(input('ratio of wave numbers k_p/k: '))
    rho_p = float(input('density of the material inside the sphere rho_p: '))
else:
    k_p_k = 1.0
    rho_p = 1.0
print()
a = float(input('radius a [mm]: '))
beta = float(input('size parameter beta: '))
kb = float(input('shift .....       kb: '))
print()
theta_p = float(input('Eulerian angle ... theta_p in [deg.]: '))
alpha = float(input('Eulerian angle ... alpha in [deg.]: '))
print()
n1cut = int(input("truncation parameter ...   n1cut: "))
print()
print()
k = beta / a 
theta_l = 180.0
theta_l1 = 91
theta = np.linspace(0.0, theta_l, theta_l1)
r_theta = theta * scp.pi / 180.

# Calculation of the T-matrix:

if para == 's':
    tm = bas.tm_s(n1cut, beta)
elif para == 'h':
    tm = bas.tm_h(n1cut, beta)
else:
    tm = bas.tm_p(n1cut, beta, a, k, k_p_k, rho_p)


# Calculation of the separation matrix:
    
sm_l_i_j = []   
for l in range(0, n1cut + 1):   
    sm_i_j = []
    for i in range(0, n1cut + 1):
        sm_j = []
        for j in range(0, n1cut + 1):
            sm_j = sm_j + [bas.SM(l, i, j, kb)]
        sm_i_j = sm_i_j + [sm_j]
    sm_l_i_j = sm_l_i_j + [sm_i_j]
sm_l_i_j = np.real(sm_l_i_j)


# Calculation of the matrix of rotation required to perform 
# the transformation into the body system:
    
dh_n_l = []
for l in range(-n1cut, n1cut + 1):   
    dh_n = []
    for n in range(0, n1cut + 1):
        dh_n = dh_n + [bas.drm(n, l, 0, 0., - theta_p, - alpha)]
    dh_n_l = dh_n_l + [dh_n]


# Calculation of the left-hand side of Eq. (2.74):

zi = 0. + 1.0j
f1 = []
for i in range(0, theta_l1):
    c_k = kb * (np.cos(theta[i] * scp.pi / 180.) * np.cos(theta_p * \
                scp.pi / 180.) + np.sin(theta[i] * scp.pi / 180.) * \
    np.sin(theta_p * scp.pi / 180.) * np.cos(alpha *scp.pi / 180.))
    c_p = kb * np.cos(theta_p * scp.pi / 180.)
    ef = np.exp(zi * (c_p - c_k))
    zw = 0.0
    for n1 in range(0, n1cut + 1):
        d_n1 = np.sqrt(4. * np.pi * (2 * n1 + 1)) * zi**n1
        Y_0_n1 = scs.sph_harm(0,n1,0.,r_theta[i])
        zw = zw + tm[n1] * d_n1 * Y_0_n1 * (-zi)**n1
    f1 = f1 + [ef * zw]
f1_r = np.real(f1)
f1_i = np.imag(f1)


# Calculation of the right-hand side of Eq. (2.74):

cmod_ls_nus = []
for ls in range(-n1cut, n1cut + 1):
    cmod_nus = []
    for nus in range(np.abs(ls), n1cut + 1):
        cmod = 0.0
        for l in range(-n1cut, n1cut + 1):
            y_l = sm_l_i_j[np.abs(l)]
            y_dh = dh_n_l[l + n1cut]
            y_l_nus = y_l[nus]
            D_rueck = bas.drm(nus, ls, l, alpha, theta_p, 0.)
            csum = 0.0
            for en in range(np.abs(l), n1cut + 1):
                pref = (-1)**(en + nus)
                d_en = np.sqrt(4. * np.pi * (2 * en + 1)) * zi**en
                D_hin = y_dh[en]
                for nu in range(np.abs(l), n1cut + 1):
                    S_rueck = y_l_nus[nu]
                    T_nu = tm[nu]
                    y_l_nu = y_l[nu]
                    S_hin = y_l_nu[en]
                    csum = csum + D_rueck * S_rueck * T_nu * S_hin * D_hin * \
                    d_en * pref
            cmod = cmod + csum        
        cmod_nus = cmod_nus + [cmod]
    cmod_ls_nus = cmod_ls_nus + [cmod_nus]

f2 = []
print()
for i in range(0, theta_l1):
    psi = 0.0
    for ls in range(-n1cut, n1cut + 1):
        c_ls = cmod_ls_nus[ls + n1cut]
        for nus in range(np.abs(ls), n1cut + 1):
            nusz = nus - np.abs(ls)
            c_ls_nus = c_ls[nusz]
            Y_ls_nus = scs.sph_harm(ls,nus,0.,r_theta[i])
            psi = psi + c_ls_nus * Y_ls_nus * (-zi)**nus
    f2 = f2 + [psi]
f2_r = np.real(f2)
f2_i = np.imag(f2)
print()


# Plot of the results:
    
plt.figure(1)
plt.subplot(121)
#plt.yscale('log') # bei bestimmten Werten besser ausschalten!
plt.plot(theta, f1_r, color = 'black', linewidth=2.0, label = "lhs")
plt.plot(theta[0], f2_r[0], 'ro', linewidth=2.0, color = 'black', label = "rhs")
for i in range(1, 9):
    plt.plot(theta[10 * i + 1], f2_r[10 * i + 1], 'ro', \
             linewidth=2.0, color = 'black')
plt.plot(theta[90], f2_r[90], 'ro', linewidth=2.0, color = 'black')
#plt.plot(theta, f2_r, color = 'r', linestyle = 'dashed', linewidth=2.0, \
#         label = "sum_phase")
plt.xlabel("scattering angle [deg]")
plt.ylabel("real")
#plt.legend()
plt.legend(loc = 'lower center')
#plt.show()


plt.subplot(122)
#plt.yscale('log') # bei bestimmten Werten besser ausschalten!
plt.plot(theta, f1_i, color = 'black', linewidth=2.0, label = "lhs")
plt.plot(theta[0], f2_i[0], 'ro', linewidth=2.0, color = 'black', label = "rhs")
for i in range(1, 9):
    plt.plot(theta[10 * i + 1], f2_i[10 * i + 1], 'ro', \
             linewidth=2.0, color = 'black')
plt.plot(theta[90], f2_i[90], 'ro', linewidth=2.0, color = 'black')
#plt.plot(theta, f2_i, color = 'r', linestyle = 'dashed', linewidth=2.0, \
#         label = "sum_phase")
plt.xlabel("scattering angle [deg.]")
plt.ylabel("imag")
#plt.legend()
plt.legend(loc = 'lower center')
plt.show()
