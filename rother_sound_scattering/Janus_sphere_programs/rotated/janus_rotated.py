#            Program "janus_rotated.py"
#
# This program calculates plane wave scattering behaviour of an
# arbitrarily rotated Janus sphere that combines the boundary 
# conditions of a sound hard, sound soft, and sound penetrable sphere 
# in three different ways.
#

print()
print()

import scipy as scp
import scipy.special as scs
import matplotlib.pyplot as plt
import basics as bas
import numpy as np
import os


# Reading scattering configuration from input file:

fobj = open("input_data_janus_rotated.txt", "r")

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
theta_j = float(z[5])
alpha = float(z[6])
theta_p = float(z[7])
ncut = int(z[8])
lcut = int(z[9])
plot_para = z[10]

k = beta / a 
kappa = rho_p / k_p_k
beta_p = k * k_p_k * a
zi = 0. + 1.0j
pref = - zi / k


# Calculation of the matrix of rotation for the transformation 
# into the body frame:
    
dh_n_l = []
for l in range(-lcut, lcut + 1):   
    dh_n = []
    for n in range(0, ncut + 1):
        dh_n = dh_n + [bas.drm(n, l, 0, 0., -theta_p, -alpha)]
    dh_n_l = dh_n_l + [dh_n]
dh_n_l = np.array(dh_n_l)


# Calculation of the expansion coefficients of the incident plane
# wave and the scattred field in the local system of the rotated 
# Janus sphere:

d_l_n = []
for l in range(-lcut, lcut + 1):
    d_n = []
    for n in range(np.abs(l), ncut + 1):
        d_inc = np.sqrt(4. * np.pi * (2 * n + 1)) * zi**n
        d_inc = d_inc * dh_n_l[l + lcut, n]
        d_n = d_n + [d_inc]
    d_l_n = d_l_n + [d_n]

c_l_n = []
for l in range(-lcut, lcut + 1):
    print(' calculating T-matrix of l-mode l = ', l)
    d_l = d_l_n[l + lcut]
    d_l = np.array(d_l)
    if para == 'hp':
        tm_l = bas.t_hp(l, beta, beta_p, theta_j, k_p_k, rho_p, ncut)
    elif para == 'sp':
        tm_l = bas.t_sp(l, beta, beta_p, theta_j, k_p_k, rho_p, ncut)
    else:
        tm_l = bas.t_hs(l, beta, theta_j, ncut)
    cv = tm_l.dot(d_l)
    c_l_n = c_l_n + [cv]


# Calculation of the differential scattering cross-section "dscross" 
# as a function of the scattering angle "theta" in the interval [0, \pi] 
# in steps of 1.0 degree:
    
theta_l = 360.0
theta_l1 = 361
theta = np.linspace(0.0, theta_l, theta_l1)
r_theta = theta * scp.pi / 180.
dscross = []
psi_s = []
for i in range(0, theta_l1):
    print(' scattering angle theta = ', theta[i])
    psi = 0.0    
    for l in range(-lcut, lcut + 1):
        cr_l = c_l_n[l + lcut]
        for n in range(np.abs(l), ncut + 1):
            nz = n - np.abs(l)
            cr = (-zi)**n * cr_l[0, nz]
            for ls in range(-n, n + 1):
                D_rueck = bas.drm(n, ls, l, alpha, theta_p, 0.)
                if i <= 180:
                    Y_ls_n = scs.sph_harm(ls,n,0.,r_theta[i])
                else:
                    Y_ls_n = scs.sph_harm(ls,n,np.pi,r_theta[360 - i])
                psi = psi + (cr * D_rueck * Y_ls_n)
    psi = pref * psi
    psi_s = psi_s + [psi]
    dscross = dscross + [psi * np.conj(psi)]                             
dscross = scp.real(dscross)


# Generation of the result file:

os.system("del dscross_janus_rotated.txt")
os.system("type nul > dscross_janus_rotated.txt")
fobj1 = open("dscross_janus_rotated.txt", "w")
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

w = scp.imag(psi_s[0])
scat_tot = 4 * scp.pi * w / k
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
