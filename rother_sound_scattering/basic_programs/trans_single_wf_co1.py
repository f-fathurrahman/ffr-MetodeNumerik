#               Program "trans_single_wf_co1.py"
#
# Transformation of the regular- and radiating wave functions "psi_{l=0,n}" 
# and "phi_{l=0,n}" for a dimensionless shift "k0b" along the positive z-axis 
# of the laboratory frame. The original and shifted wave functions are 
# calculated at different angles "theta1" in [0, \pi] in steps of 
# "180/w_num_m" degrees, and for a fixed but dimensionless distance "k_0r_1" in 
# the shifted co.-system!
#


print()
print()
print(" --- Translation of a regular- and radiating  wave function ...")
print("       ... along the positive z-axis. k0b >= 2 * k0r1 ! --- ")
print()
print()

import numpy as np
import scipy as scp
import scipy.special as scs
import basics as bas
import matplotlib.pyplot as plt

    
# Input order of Legendre polynomial, the transformation parameters, 
# and the truncation paerameter:
    
en = int(input('order of Legendre polynomial  n: '))
k0r1 = float(input('dimensionless distance in the shifted system ... k0r1: '))
a1 = 2. * k0r1
a1 = np.str(a1)
a0 = "dimensionless shift along positive z-axis (k0b >= "
a2 = "!)  k0b: "
a = a0 + a1 + a2
k0b = float(input(a))
nucut = int(input('truncation parameter ... nucut: '))
print()
    
w_num = 37         #  number of angles "\theta1" in [0°,180°]!
w_num_m = w_num - 1
zi = 0. + 1.0j
theta1_o1 = np.linspace(0.0, 180.0, w_num)
ctheta1_o1 = np.cos(theta1_o1 * scp.pi / 180.)
xa = np.arange(w_num)
xachse = np.arange(w_num)


# calculation of the separation matrix:

y_os = []
for s_sum in range(0, nucut + 1):
    y_os = y_os + [bas.SM(0, en, s_sum, k0b)] 
y_o1 = np.real(y_os)


# calculation of the dimensionless distances "kr" and angles "theta" in
# the laboratory frame:
    
theta1 = 180. - theta1_o1
ctheta1 = np.cos(theta1 * scp.pi / 180.)
theta = []
ctheta = []
k0r = []    
for i in range(0,w_num):
    k0r = k0r + [np.sqrt(k0b**2 + k0r1**2 - 2. * k0b * k0r1 * ctheta1[i])]
    if i == 0 or i == w_num_m:
        arg_cos = 1.0
    else:
        arg_cos = -(k0r1**2 - k0r[i]**2 - k0b**2) / (2. * k0r[i] * k0b)
    theta = theta + [180.0 * np.arccos(arg_cos) / scp.pi]
    ctheta = ctheta + [np.cos(theta[i] * scp.pi / 180.)]
ctheta = np.array(ctheta)
theta = np.array(theta)             
kr = np.array(k0r)


# calculation of the wave functions in the laboratory frame:

psi_n = []
psi_n_trans = []
phi_n = []
phi_n_trans = []
for j in range(0, w_num):
    pn0 = scs.lpn(en,ctheta[j])
    pn1 = pn0[0]
    ysh = np.sqrt((2. * en + 1) / 4. / np.pi) * float(pn1[en])
    u = scs.spherical_jn(en,kr[j])
    v0 = scs.spherical_yn(en,kr[j])
    v = u + zi * v0
    psi_n = psi_n + [u * ysh]
    phi_n = phi_n + [v * ysh] 

     
# transformation of the wave functions of the laboratory frame:  
    
    pn0_o1 = scs.lpn(nucut,ctheta1_o1[j])
    pn1_o1 = pn0_o1[0]
    tpsi = 0.
    tphi = 0.
    for nu in range(0, nucut + 1):
        ysh_o1 = np.sqrt((2. * nu + 1) / 4. / np.pi) * pn1_o1[nu]
        u_o1 = scs.spherical_jn(nu,k0r1)
        v0_o1 = scs.spherical_yn(nu,k0r1)
        v_o1 = u_o1 + zi * v0_o1
        psi_nu_xr1 = u_o1 * ysh_o1
        phi_nu_xr1 = v_o1 * ysh_o1
        tpsi = tpsi + y_o1[nu] * psi_nu_xr1
        tphi = tphi + y_os[nu] * psi_nu_xr1
    tpsi = (-1)**en * tpsi
    tphi = (-1)**en * tphi
    psi_n_trans = psi_n_trans + [tpsi]
    phi_n_trans = phi_n_trans + [tphi]

psi_n = np.array(psi_n)
psi_n_trans = np.array(psi_n_trans)  
phi_n = np.array(phi_n)
phi_n_trans = np.array(phi_n_trans)    


# Plot of the results:
    
plt.plot(1)
plt.plot(theta1_o1, psi_n, 'o', color = 'black', linewidth = 2.0, label = \
         "\psi_n")
plt.plot(theta1_o1, psi_n_trans, color = 'black', linewidth = 2.0, label = \
         "\psi_n_trans")
plt.legend()
plt.xlabel("$\\theta_1$ (degree)", fontsize=16)
plt.show()

plt.plot(2)
plt.plot(theta1_o1, np.real(phi_n), 'o', color = 'black', linewidth = 2.0, label = \
         "real(\phi_n)")
plt.plot(theta1_o1, np.real(phi_n_trans), color = 'black', linewidth = 2.0, \
         label = "real(\phi_n_trans)")
plt.legend()
plt.xlabel("$\\theta_1$ (degree)", fontsize=16)
plt.show()

plt.plot(3)
plt.plot(theta1_o1, np.imag(phi_n), 'o', color = 'black', linewidth = 2.0, label = \
         "imag(\phi_n)")
plt.plot(theta1_o1, np.imag(phi_n_trans), color = 'black', linewidth = 2.0, \
         label = "imag(\phi_n_trans)")
plt.legend()
plt.xlabel("$\\theta_1$ (degree)", fontsize=16)
plt.show()
