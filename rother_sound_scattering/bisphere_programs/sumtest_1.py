#            Program "sumtest_1.py"
#
# This program provides a numerical test of Eqs. (4.31) and (4.32) in Chpt. 4.
#

print()
print()
print("                      --- Sum test 1! ---")
print()
print()

import numpy as np
import scipy as scp
import scipy.special as scs
import matplotlib.pyplot as plt
import basics as bas


# Input:
    
kb = float(input(' real; shift ... k_0*b : '))
m1 = int(input(' integer; mode ... m0 : '))
l1 = int(input(' integer; mode ... l1 : '))
ncut = int(input(' integer; truncation parameter ... ncut : '))


# Calculation of the analytical function and the approximation:

zi = 0. + 1.0j 
y_m0_m1 = []   
for i in range(0, ncut + 1):   
    y_os = []
    for j in range(0, ncut + 1):
        y_os = y_os + [bas.SM(l1, i, j, kb)]
    y_m0_m1 = y_m0_m1 + [y_os]
y_m0_m1 = np.real(y_m0_m1)
S_m1 = y_m0_m1[m1]

theta = np.linspace(0.0, 180.0, 181)
r_theta = theta * scp.pi / 180.

sum_phase = []
ana_phase = []
for i in range(0, 181):
    c_k = kb * np.cos(r_theta[i])
    Y_l1_m1 = scs.sph_harm(l1,m1,0.,r_theta[i])
    sum = 0.0
    abs_l1 = np.abs(l1)
    for m1s in range(abs_l1, ncut + 1):
        Y_l1_m1s = scs.sph_harm(l1,m1s,0.,r_theta[i])
        sum = sum + zi**(m1s + m1) * S_m1[m1s] * Y_l1_m1s
    sum_phase = sum_phase + [sum]
    ana = np.exp(- zi * c_k) * Y_l1_m1
    ana_phase = ana_phase + [ana]

f1_r = np.real(ana_phase)
f2_r = np.real(sum_phase)

f1_i = np.imag(ana_phase)
f2_i = np.imag(sum_phase)

    
plt.figure(1)
plt.subplot(121)
#plt.yscale('log')
plt.plot(theta, f1_r, color = 'black', linewidth=2.0, label = "analytical")
plt.plot(theta[0], f2_r[0], 'ro', linewidth=2.0, color = 'black', label = "sum")
for i in range(1, 18):
    plt.plot(theta[10 * i + 1], f2_r[10 * i + 1], 'ro', \
             linewidth=2.0, color = 'black')
plt.plot(theta[180], f2_r[180], 'ro', linewidth=2.0, color = 'black')
#plt.plot(theta, f2_r, color = 'black', linestyle = 'dashed', linewidth=2.0, \
#         label = "sum")
plt.xlabel("scattering angle [deg.]", fontsize=16)
plt.ylabel("real", fontsize=16)
plt.legend(loc = 'lower center')
 
plt.subplot(122)
#plt.yscale('log')
plt.plot(theta, f1_i, color = 'black', linewidth=2.0, label = "analytical")
plt.plot(theta[0], f2_i[0], 'ro', linewidth=2.0, color = 'black', label = "sum")
for i in range(1, 18):
    plt.plot(theta[10 * i + 1], f2_i[10 * i + 1], 'ro', \
             linewidth=2.0, color = 'black')
plt.plot(theta[180], f2_i[180], 'ro', linewidth=2.0, color = 'black')
#plt.plot(theta, f2_i, color = 'black', linestyle = 'dashed', linewidth=2.0, \
#         label = "sum")
plt.xlabel("scattering angle [deg.]", fontsize=16)
plt.ylabel("imag", fontsize=16)
plt.legend(loc = 'upper center')
plt.show()
