#            Program "shifted_sphere.py"
#
# Acoustic plane wave scattering on an arbitrarily shifted sphere. 
# Calculation is not performed by use of a translation of the laboratory 
# frame into the body frame but by using only a rotation and a shift 
# along the new z-axis after the rotation! Note, that, if a penetrable 
# sphere is considered, the density of the material "rho_0" outside the 
# spheres is always given by 1.

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

fobj = open("input_data_shifted_sphere.txt", "r")

z = []
for line in fobj:
    line = line.strip()
    arr = line.split("= ")
    wert = str(arr[1])
    z = z + [wert]
fobj.close()

para = z[0]
k_p_k = float(z[1])
rho_p = float(z[2])
a = float(z[3])
beta = float(z[4])
alpha = float(z[5])
theta_p = float(z[6])
b = float(z[7])
ncut = int(z[8])
plot_para = z[9]

nangle = 91  # defines resolution in the scattering plane [0, \pi]
k0 = beta / a
k0b = k0 * b
zi = 0. + 1.0j


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
sm_l_i_j = np.real(sm_l_i_j)


# Calculation of the matrix of rotation required to perform 
# the transformation into the body system:
    
dh_n_l = []
for l in range(-ncut, ncut + 1):   
    dh_n = []
    for n in range(0, ncut + 1):
        dh_n = dh_n + [bas.drm(n, l, 0, 0., - theta_p, - alpha)]
    dh_n_l = dh_n_l + [dh_n]


# Calculation of the T-matrix:

if para == 's':
    tm = bas.tm_s(ncut, beta)
elif para == 'h':
    tm = bas.tm_h(ncut, beta)
else:
    tm = bas.tm_p(ncut, beta, a, k0, k_p_k, rho_p)


# Calculation of the differential scattering cross-section "dscross" 
# in the laboratory frame as a function of the scattering angle "theta" 
# in the interval [0, \pi] in steps depending on "nangle":
    
cmod_ls_nus = []
for ls in range(-ncut, ncut + 1):
    cmod_nus = []
    for nus in range(np.abs(ls), ncut + 1):
        cmod = 0.0
        for l in range(-ncut, ncut + 1):
            y_l = sm_l_i_j[np.abs(l)]
            y_dh = dh_n_l[l + ncut]
            y_l_nus = y_l[nus]
            D_rueck = bas.drm(nus, ls, l, alpha, theta_p, 0.)
            csum = 0.0
            for en in range(np.abs(l), ncut + 1):
                pref = (-1)**en * zi**nus
                d_en = np.sqrt(4. * np.pi * (2 * en + 1)) * zi**en
                D_hin = y_dh[en]
                for nu in range(np.abs(l), ncut + 1):
                    S_rueck = y_l_nus[nu]
                    T_nu = tm[nu]
                    y_l_nu = y_l[nu]
                    S_hin = y_l_nu[en]
                    csum = csum + D_rueck * S_rueck * T_nu * S_hin * D_hin * \
                    d_en * pref
            cmod = cmod + csum        
        cmod_nus = cmod_nus + [cmod]
    cmod_ls_nus = cmod_ls_nus + [cmod_nus]

theta = np.linspace(0.0, 180.0, nangle)
r_theta = theta * scp.pi / 180.
pf = - zi / k0
dscross = []
psi_s = []
print()
for i in range(0, nangle):
    psi = 0.0
    for ls in range(-ncut, ncut + 1):
        c_ls = cmod_ls_nus[ls + ncut]
        for nus in range(np.abs(ls), ncut + 1):
            nusz = nus - np.abs(ls)
            c_ls_nus = c_ls[nusz]
            Y_ls_nus = scs.sph_harm(ls,nus,0.,r_theta[i])
            psi = psi + c_ls_nus * Y_ls_nus
    psi = pf * psi
    psi_s = psi_s + [psi]
    dscross = dscross + [psi * np.conj(psi)]
dscross = np.real(dscross)
print()


# Generation of the result file:

os.system("del dscross_shifted_sphere.txt")
os.system("type nul > dscross_shifted_sphere.txt")
fobj1 = open("dscross_shifted_sphere.txt", "w")
nl = "\n"
for i in range(0, nangle):
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
scat_tot = 4 * scp.pi * w / k0
print()
print("total scattering cross-section: scat_tot = ", scat_tot)


# Plot of "dscross":

if plot_para == 'lg':
    plt.yscale('log')
    plt.plot(theta, dscross, color='black', linewidth=2.0)
else:    
    plt.plot(theta, dscross, color='black', linewidth=2.0)
plt.xlabel("scattering angle [degree]", fontsize=16)
plt.ylabel("diff. scat. cross-sect.", fontsize=16)
plt.show()
