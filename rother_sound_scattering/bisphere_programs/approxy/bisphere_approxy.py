#            Program "bisphere_approxy.py"
#
# This program approximates the plane wave scattering behaviour of any 
# combination of two acoustically soft, hard or penetrable spheres that
# are arbitrarily oriented in the laboratory frame! This approximation 
# neglects any interaction between the spheres and consideres only the 
# interference term. The configuration is fixed by the distance "b", and 
# the Eulerian angles "alpha" and "\theta_p". 
#

print()
print()

import numpy as np
import scipy as scp
import scipy.special as scs
import matplotlib.pyplot as plt
import basics as bas
import os


# Read scattering configuration from input file:

fobj = open("input_data_bisphere_approxy.txt", "r")

z = []
for line in fobj:
    line = line.strip()
    arr = line.split("= ")
    wert = str(arr[1])
    z = z + [wert]
fobj.close()

para0 = z[0]
para1 = z[1]
k_p_k_0 = float(z[2])
rho_p_0 = float(z[3])
k_p_k_1 = float(z[4])
rho_p_1 = float(z[5])
a0 = float(z[6])
a1 = float(z[7])
beta_0 = float(z[8])
b = float(z[9])
alpha = float(z[10])
theta_p = float(z[11])
ppm = int(z[12])
plot_para = z[13]


# Fixing angular resolution and truncation parameter:

if ppm == 1:
    theta_l = 360.0
    theta_l1 = 361
    theta = np.linspace(0.0, theta_l, theta_l1)
    ctheta = np.cos(theta * scp.pi / 180.)
    ctheta[0] = 1.0
    ctheta[180] = -1.0
else:
    theta_l = 180.0
    theta_l1 = 181
    theta = np.linspace(0.0, theta_l, theta_l1)
    ctheta = np.cos(theta * scp.pi / 180.)
    ctheta[0] = 1.0

n0cut = int(beta_0 + 9.) # may be changed for larger size parameters!
k = beta_0 / a0
kb = k * b
beta_1 = k * a1
n1cut = int(beta_1 + 9.) # may be changed for larger size parameters!


# Calculation of the T-matrices:

if (para0 == 's' and para1 == 's'):
    tm_0 = bas.tm_s(n0cut, beta_0)
    tm_1 = bas.tm_s(n1cut, beta_1)
elif (para0 == 's' and para1 == 'h'):
    tm_0 = bas.tm_s(n0cut, beta_0)
    tm_1 = bas.tm_h(n1cut, beta_1)
elif (para0 == 's' and para1 == 'p'):
    tm_0 = bas.tm_s(n0cut, beta_0)
    tm_1 = bas.tm_p(n1cut, beta_1, a1, k, k_p_k_1, rho_p_1)
elif (para0 == 'h' and para1 == 's'):
    tm_0 = bas.tm_h(n0cut, beta_0)
    tm_1 = bas.tm_s(n1cut, beta_1)
elif (para0 == 'h' and para1 == 'h'):
    tm_0 = bas.tm_h(n0cut, beta_0)
    tm_1 = bas.tm_h(n1cut, beta_1)
elif (para0 == 'h' and para1 == 'p'):
    tm_0 = bas.tm_h(n0cut, beta_0)
    tm_1 = bas.tm_p(n1cut, beta_1, a1, k, k_p_k_1, rho_p_1)
elif (para0 == 'p' and para1 == 's'):
    tm_0 = bas.tm_p(n0cut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_s(n1cut, beta_1)
elif (para0 == 'p' and para1 == 'h'):
    tm_0 = bas.tm_p(n0cut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_h(n1cut, beta_1)
else:
    tm_0 = bas.tm_p(n0cut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_p(n1cut, beta_1, a1, k, k_p_k_1, rho_p_1)

    
# Calculation of the differential scattering cross-section as a 
# function of the scattering angle "theta" in the interval [0, \pi] 
# or [0, 2 * \pi] in steps of 1 degree:

zi = 0. + 1.0j
pref = -zi / k
dscross = []
psi_s = []
for i in range(0, theta_l1):
    psi0_1 = 0.0
    psi1_1 = 0.0
    y0 = scs.lpn(n0cut,ctheta[i])
    y1 = scs.lpn(n1cut,ctheta[i])
    y0_1 = y0[0]
    y1_1 = y1[0]
    for n0 in range(0, n0cut + 1):
        d_n0 = (2 * n0 + 1)
        sum_n0 = d_n0 * y0_1[n0] * tm_0[n0]
        psi0_1 = psi0_1 + sum_n0    
    for n1 in range(0, n1cut + 1):
        d_n1 = (2 * n1 + 1)
        sum_n1 = d_n1 * y1_1[n1] * tm_1[n1]
        psi1_1 = psi1_1 + sum_n1
    psi0_1 = pref * psi0_1
    psi1_1 = pref * psi1_1
    c_k = kb * (np.cos(theta[i] * scp.pi / 180.) * np.cos(theta_p * \
                scp.pi / 180.) + np.sin(theta[i] * scp.pi / 180.) * \
    np.sin(theta_p * scp.pi / 180.) * np.cos(alpha *scp.pi / 180.))               
    c_p = kb * np.cos(theta_p * scp.pi / 180.)
    psi1_1 = psi1_1 * np.exp(zi * (c_p - c_k))
    psi = psi0_1 + psi1_1
    psi_s = psi_s + [psi]
    dscross = dscross + [psi * np.conj(psi)]
dscross = np.real(dscross)


# Generation of the result file:

os.system("del dscross_bisphere_approxy.txt")
os.system("type nul > dscross_bisphere_approxy.txt")
fobj1 = open("dscross_bisphere_approxy.txt", "w")
nl = "\n"
for i in range(0, theta_l1):
    a1 = str(theta[i])
    a2 = " = "
    a3 = str(dscross[i]) + nl
    a = a1 + a2 + a3
    fobj1.write(a)
fobj1.close()


# Calculation of the total scattering cross-section "scat_tot" by use 
# of the optical theorem:

print()
print()
print('Results: ')
w = np.imag(psi_s[0])
scat_tot = 4 * scp.pi * w / k
print()
print("total scattering cross-section: scat_tot = ", scat_tot)
print()


# Plot of "dscross":

if plot_para == 'lg':
    plt.yscale('log')
    plt.plot(theta, dscross, color = 'black', linewidth=2.0)
else:    
    plt.plot(theta, dscross, color = 'black', linewidth=2.0)
plt.xlabel("scattering angle [degree]",fontsize=16)
plt.ylabel("diff. scat. cross-sect.",fontsize=16)
plt.show()
