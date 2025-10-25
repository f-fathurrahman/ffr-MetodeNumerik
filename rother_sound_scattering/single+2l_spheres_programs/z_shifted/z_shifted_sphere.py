#            Program "z_shifted_sphere.py"
#
# This program calculates the plane wave scattering behaviour of an 
# acoustically soft, hard or penetrable sphere (homogeneous Dirichlet-, 
# von Neumann- or mixed condition) shifted along the z-axis of the 
# laboratory frame. Note, that the density of the material "rho_0" outside 
# the spheres is always given by 1. This program uses the expansion of 
# the primary incident plane wave in the laboratory frame and its 
# transformation into the shifted system afterwards.
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

fobj = open("input_data_z_shifted_sphere.txt", "r")

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
b = float(z[5])
ncut = int(z[6])
plot_para = z[7]

k0 = beta / a
zi = 0. + 1.0j
k0b = k0 * b


# Calculation of the T-matrix:

if para == 's':
    tm = bas.tm_s(ncut, beta)
elif para == 'h':
    tm = bas.tm_h(ncut, beta)
else:
    tm = bas.tm_p(ncut, beta, a, k0, k_p_k, rho_p)


# Calculation of the separation matrix:
    
y_os_n_nu = []   
for i in range(0, ncut + 1):   
    y_os_nu = []
    for j in range(0, ncut + 1):
        y_os_nu = y_os_nu + [bas.SM(0, i, j, k0b)]
    y_os_n_nu = y_os_n_nu + [y_os_nu]
y_o1_n_nu = np.real(y_os_n_nu)


# Calculation of the expansion coefficients in the laboratory frame:
    
def c_0_n_fun(nu_1, beta, kb, ncut):
    c = 0.0
    y_o_n = y_o1_n_nu[nu_1]
    for nu in range(0, ncut + 1):
        pref = (-1)**nu_1
        d_nu = np.sqrt(4. * np.pi * (2 * nu + 1)) * zi**nu
        S_rueck = y_o_n[nu]
        c = c + pref * S_rueck * tm[nu] * d_nu
    c = c * np.exp(zi * k0b)
    return c


# Calculation of the differential scattering cross-section "dscross" 
# as a function of the scattering angle "theta" in the interval [0, \pi] 
# in steps of 1 degree:
    
theta = np.linspace(0.0, 180.0, 181)
r_theta = theta * scp.pi / 180.
pf = - zi / k0
dscross = []
psi_s = []
print()
for i in range(0, 181):
    psi = 0.0
    for nu_1 in range(0, ncut + 1):
        c_ns = c_0_n_fun(nu_1, beta, k0b, ncut)
        Y_0_nu1 = scs.sph_harm(0,nu_1,0.,r_theta[i])
        sum = (-zi)**nu_1 * c_ns * Y_0_nu1
        psi = psi + sum
    psi = pf * psi
    psi_s = psi_s + [psi]
    dscross = dscross + [psi * np.conj(psi)]
dscross = np.real(dscross)
print()


# Generation of the result file:

os.system("del dscross_z_shifted_sphere.txt")
os.system("type nul > dscross_z_shifted_sphere.txt")
fobj1 = open("dscross_z_shifted_sphere.txt", "w")
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



