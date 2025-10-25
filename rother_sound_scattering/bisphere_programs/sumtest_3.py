#            Program "sumtest_3.py"
#
# This program provides a test of Eq. (2.75) of the manuscript.
#

print()
print()
print("         Test of Eq. (2.75) of the manuscript!")
print()
print()

import numpy as np
import scipy as scp
import scipy.special as scs
import matplotlib.pyplot as plt
import basics as bas
    

# Input parameter:
    
n1 = int(input("order .................      n1: "))
kb = float(input('shift ..............         kb: '))
n1cut = int(input("truncation parameter ...     n1cut: "))
print()

theta_p = 90.
theta_l = 180.0
theta_l1 = 91
theta = np.linspace(0.0, theta_l, theta_l1)
r_theta = theta * scp.pi / 180.
zi = 0. + 1.0j
d_n1 = np.sqrt(4. * np.pi * (2 * n1 + 1)) * zi**n1


# Calculation of the left-hand side of Eq. (2.75):

f1 = []
for i in range(0, theta_l1):
    c_k = kb * np.cos(np.abs(theta[i] - theta_p) * scp.pi / 180.)
    c_p = kb * np.cos(theta_p * scp.pi / 180.)
    ef = np.exp(zi * c_k)
    Y_0_n1 = scs.sph_harm(0,n1,0.,r_theta[i])
    zw = ef * d_n1 * Y_0_n1
    f1 = f1 + [zw]
f1_r = np.real(f1)
f1_i = np.imag(f1)


# Calculation of the right-hand side of Eq. (2.75):

f2_tl1_nu1 = []
for tl1 in range(-n1cut, n1cut + 1):
    f2_nu1 = []
    for nu1 in range(np.abs(tl1), n1cut + 1):
        f2_zs = 0.0
        for l1 in range (-n1, n1 + 1):                
            D_rueck = bas.drm(nu1, tl1, l1, 0., theta_p, 0.)
            S_rueck = bas.SM(l1, nu1, n1, kb)
            S_rueck = np.real(S_rueck)                
            for n0 in range(0, n1cut + 1):
                vf1 = zi**(nu1 + n1)
                vf2 = (-1)**n0
                d_n0 = np.sqrt(4. * np.pi * (2 * n0 + 1)) * zi**n0
                vf3 = vf1 * vf2 * d_n0
                D_hin = bas.drm(n0, l1, 0, 0., -theta_p, 0.)
                S_hin = bas.SM(l1, n0, n1, kb)
                S_hin = np.real(S_hin)
                f2_zs = f2_zs + vf3 * S_hin * D_hin * S_rueck * D_rueck        
        f2_nu1 = f2_nu1 + [f2_zs]
    f2_tl1_nu1 = f2_tl1_nu1 + [f2_nu1]
 
f2 = []
for i in range(0, theta_l1):
    fw = 0.0
    for tl1 in range(-n1cut, n1cut + 1):
        fw_zs_tl1 = f2_tl1_nu1[tl1 + n1cut]
        for nu1 in range(np.abs(tl1), n1cut + 1):
            nu1z = nu1 - np.abs(tl1)
            Y_tl1_nu1 = (-1)**tl1 * scs.sph_harm(tl1,nu1,0.,r_theta[i])
            fw_zs = fw_zs_tl1[nu1z] * Y_tl1_nu1
            fw = fw + fw_zs
    f2 = f2 + [fw]
f2_r = np.real(f2)
f2_i = np.imag(f2)


# Plot of the results

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
