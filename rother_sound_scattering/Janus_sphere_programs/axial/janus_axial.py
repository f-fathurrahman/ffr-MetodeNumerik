#            Program "janus_axial.py"
#
# This program calculates plane wave scattering behaviour of an
# axisymmetrically oriented Janus sphere that combines the boundary 
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

fobj = open("input_data_janus_axial.txt", "r")

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
ncut = int(z[6])
plot_para = z[7]

k = beta / a 
kappa = rho_p / k_p_k
beta_p = k * k_p_k * a
theta = np.linspace(0.0, 180.0, 181)
r_theta = theta * scp.pi / 180.
zi = 0. + 1.0j
pref = - zi / k

d = []
for n in range(0, ncut +1):
    d = d + [np.sqrt(4. * np.pi * (2 * n + 1)) * zi**n]


# Calculation of the T-matrix:

if para == 'hp':
    tm = bas.t_hp(0, beta, beta_p, theta_j, k_p_k, rho_p, ncut)
elif para == 'sp':
    tm = bas.t_sp(0, beta, beta_p, theta_j, k_p_k, rho_p, ncut)
else:
    tm = bas.t_hs(0, beta, theta_j, ncut)


# Calculation of the scattering coefficients:
    
c = tm.dot(d)   


# Calculation of the differential scattering cross-section "dscross" 
# as a function of the scattering angle "theta" in the interval [0, \pi] 
# in steps of 1.0 degree:
    
dscross = []
psi_s = []
for i in range(0, 181):
    psi = 0.0
    for n in range(0, ncut +1):
        Y_0_n = scs.sph_harm(0,n,0.,r_theta[i])
        sum_n = c[0,n] * (-zi)**n * Y_0_n
        psi = psi + sum_n
    psi = pref * psi
    psi_s = psi_s + [psi]
    dscross = dscross + [psi * np.conj(psi)]
dscross = np.real(dscross)


# Generation of the result file:

os.system("del dscross_janus_axial.txt")
os.system("type nul > dscross_janus_axial.txt")
fobj1 = open("dscross_janus_axial.txt", "w")
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

w = np.imag(psi_s[0])
scat_tot = 4 * scp.pi * w / k
print()
print("total scattering cross-section: scat_tot = ", scat_tot)


# Plot of "dscross":

if plot_para == 'lg':
    plt.yscale('log')
    plt.plot(theta, dscross, color='black', linewidth=2.0)
else:    
    plt.plot(theta, dscross, color='black', linewidth=2.0)
plt.xlabel("scattering angle [degree]", fontsize = 16)
plt.ylabel("diff. scat. cross-sect.", fontsize = 16)
plt.show()
