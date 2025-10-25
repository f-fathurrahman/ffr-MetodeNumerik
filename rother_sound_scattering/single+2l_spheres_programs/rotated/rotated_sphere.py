#            Program "rotated_sphere.py"
#
# This program calculates the plane wave scattering behaviour of an 
# acoustically soft, hard or penetrable sphere (homogeneous Dirichlet-, 
# von Neumann- or mixed condition) arbitrarily rotated in the 
# laboratory frame. Note that, if a penetrable sphere is considered, 
# the density of the material "rho_0" outside the spheres is always given 
# by 1.
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

fobj = open("input_data_rotated_sphere.txt", "r")

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
ncut = int(z[7])
plot_para = z[8]

k0 = beta / a
zi = 0. + 1.0j


# Calculation of the T-matrix:

if para == 's':
    tm = bas.tm_s(ncut, beta)
elif para == 'h':
    tm = bas.tm_h(ncut, beta)
else:
    tm = bas.tm_p(ncut, beta, a, k0, k_p_k, rho_p)


# Calculation of the expansion coefficients in the laboratory frame:

c_n_l = []
for n in range(0, ncut + 1):
    d_n = np.sqrt(4. * np.pi * (2. * n + 1.)) * zi**n
    t_n = tm[n]
    c_l = []
    for l in range(-n, n + 1):
        c = 0.0
        for l1 in range(-n, n + 1):
            D_hin = bas.drm(n, l1, 0, 0., - theta_p, - alpha)
            D_rueck = bas.drm(n,l,l1, alpha, theta_p, 0.)
            c = c + d_n * t_n * D_hin * D_rueck
        c_l = c_l + [c]
    c_n_l = c_n_l + [c_l]


# Calculation of the differential scattering cross-section "dscross" 
# as a function of the scattering angle "theta" in the interval [0, \pi] 
# in steps of 1 degree:

theta = np.linspace(0.0, 180.0, 181)
r_theta = theta * scp.pi / 180.
ctheta = np.cos(r_theta)
pref = - zi / k0
dscross = []
psi_s = []
for i in range(0, 181):
    psi = 0.0
    for n in range(0, ncut + 1):
        c_n = c_n_l[n]
        for l in range(-n, n + 1):
            c_nl = c_n[l + n]
            Y_nl = scs.sph_harm(l,n,0.,r_theta[i])
            psi = psi + (-zi)**n * c_nl * Y_nl
    psi = pref * psi
    psi_s = psi_s + [psi]
    dscross = dscross + [psi * np.conj(psi)]
dscross = np.real(dscross)


# Generation of the result file:

os.system("del dscross_rotated_sphere.txt")
os.system("type nul > dscross_rotated_sphere.txt")
fobj1 = open("dscross_rotated_sphere.txt", "w")
nl = "\n"
for i in range(0, 181):
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
