#               Program "trans_single_wf_co.py"
#
# Transformation of the regular- and radiating wave functions "psi_{l=0,n}" 
# and "phi_{l=0,n}" for a dimensionless shift "k0b" along the positive z-axis 
# of the laboratory frame. The original and shifted wave functions are 
# calculated at different angles "theta" in [0, \pi] in steps of 
# "180/w_num_m" degrees for a fixed but dimensionless distance "k_0r" in 
# the laboratory system!
#

print()
print()
print(" --- Translation of a regular- and radiating  wave ...")
print("     ... function along the positive z-axis! --- ")

print()
print()

import numpy as np
import scipy as scp
import scipy.special as scs
import basics as bas
import matplotlib.pyplot as plt


# Input order of Legendre polynomial, the transformation parameters, 
# and the truncatuion parameter:
    
en = int(input('order of Legendre polynomial ... n: '))
k0r = float(input('dimensionless distance in laboratory frame ... k0r: '))
k0b = float(input('dimensionless shift along z-axis ... k0b: '))
nucut = int(input('truncation parameter ... nucut: '))

w_num = 37       #  number of angles "\theta" in [0°,180°]!
w_num_m = w_num - 1
zi = 0. + 1.0j
theta = np.linspace(0.0, 180.0, w_num)
ctheta = np.cos(theta * scp.pi / 180.)


# Calculation of the dimensionless distances "kr1" and angles "theta1" in
# the shifted c.-system:
    
theta1 = []
ctheta1 = []
k0r1 = []    
for i in range(0,w_num):
    k0r1 = k0r1 + [np.sqrt(k0b**2 + k0r**2 - 2. * k0b * k0r * ctheta[i])]
    if i == 0:
        arg_cos = -1.0
    elif i == w_num_m:
        arg_cos = 1.0
    else:
        arg_cos = -(k0r**2 - k0r1[i]**2 - k0b**2) / (2. * k0r1[i] * k0b)
    theta1 = theta1 + [180. * ( 1. - np.arccos(arg_cos) / scp.pi)]
    if k0b <= k0r:
        theta1[0] = 0.0
    else:
        theta1[0] = 180.0
    ctheta1 = ctheta1 + [np.cos(theta1[i] * scp.pi / 180.)]
ctheta1 = np.array(ctheta1)
theta1 = np.array(theta1)              
kr1 = np.array(k0r1)


# Calculation of the separation matrix:

y_os = []
for s_sum in range(0, nucut + 1):
    y_os = y_os + [bas.SM(0, en, s_sum, k0b)] 
y_o1 = np.real(y_os)


# calculation of the wave functions in the laboratory frame:

psi_n = []
psi_n_trans = []
phi_n = []
phi_n_trans = []
for j in range(0, w_num):
    pn0 = scs.lpn(en,ctheta[j])
    pn1 = pn0[0]
    ysh = np.sqrt((2. * en + 1) / 4. / np.pi) * float(pn1[en])
    u = scs.spherical_jn(en,k0r)
    v0 = scs.spherical_yn(en,k0r)
    v = u + zi * v0
    psi_n = psi_n + [u * ysh]
    phi_n = phi_n + [v * ysh] 

     
# transformation of the wave functions of the laboratory frame: 
     
    pn0_o1 = scs.lpn(nucut,ctheta1[j])
    pn1_o1 = pn0_o1[0]
    tpsi = 0.
    tphi = 0.
    for nu in range(0, nucut + 1):
        ysh_o1 = np.sqrt((2. * nu + 1) / 4. / np.pi) * float(pn1_o1[nu])
        u_o1 = scs.spherical_jn(nu,kr1[j])
        v0_o1 = scs.spherical_yn(nu,kr1[j])
        v_o1 = u_o1 + zi * v0_o1
        psi_nu_xr1 = u_o1 * ysh_o1
        phi_nu_xr1 = v_o1 * ysh_o1
        tpsi = tpsi + y_o1[nu] * psi_nu_xr1
        if kr1[j] < k0b:
            tphi = tphi + y_os[nu] * psi_nu_xr1
        else:
            tphi = tphi + y_o1[nu] * phi_nu_xr1
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
plt.plot(theta, psi_n, 'o', color = 'black', linewidth = 2.0, label = \
         "\psi_n")
plt.plot(theta, psi_n_trans, color = 'black', linewidth = 2.0, label = \
         "\psi_n_trans")
plt.legend()
plt.xlabel("$\\theta$ (degree)", fontsize=16)
plt.show()

plt.plot(2)
plt.plot(theta, np.real(phi_n), 'o', color = 'black', linewidth = 2.0, label = \
         "real(\phi_n)")
plt.plot(theta, np.real(phi_n_trans), color = 'black', linewidth = 2.0, \
         label = "real(\phi_n_trans)")
plt.legend()
plt.xlabel("$\\theta$ (degree)", fontsize=16)
plt.show()

plt.plot(3)
plt.plot(theta, np.imag(phi_n), 'o', color = 'black', linewidth = 2.0, label = \
         "imag(\phi_n)")
plt.plot(theta, np.imag(phi_n_trans), color = 'black', linewidth = 2.0, \
         label = "imag(\phi_n_trans)")
# plt.legend()
plt.legend(loc = 'lower right')
plt.xlabel("$\\theta$ (degree)", fontsize=16)
plt.show()
