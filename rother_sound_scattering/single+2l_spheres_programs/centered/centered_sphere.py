#            Program "centered_sphere.py"
#
# This program calculates the plane wave scattering behaviour of an 
# acoustically soft, hard, or penetrable sphere (homogeneous Dirichlet-, 
# von Neumann- or mixed condition) centered in the laboratory frame. 
# Note that the density of the material "rho_0" outside the spheres is 
# always given by 1, and that the truncation parameter "ncut" is fixed 
# to "ncut = beta + 5".
#

print()
print()

import scipy as scp
import scipy.special as scs
import matplotlib.pyplot as plt
import basics as bas
import os


# Read scattering configuration from input file:

fobj = open("input_data_centered_sphere.txt", "r")

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
plot_para = z[5]


# Fixing the truncation parameter "ncut" and the wave number "k0":
    
ncut = int(beta + 5.) # must possibly be increased for higher size parameters
k0 = beta / a 


# Calculation of the T-matrix:

if para == 's':
    tm = bas.tm_s(ncut, beta)
elif para == 'h':
    tm = bas.tm_h(ncut, beta)
else:
    tm = bas.tm_p(ncut, beta, a, k0, k_p_k, rho_p)


# Calculation of the differential scattering cross-section "dscross" 
# as a function of the scattering angle "theta" in the interval [0, \pi] 
# in steps of 0.5 degree:
    
theta = scp.linspace(0.0, 180.0, 181)
ctheta = scp.cos(theta * scp.pi / 180.)
zi = 0. + 1.0j
pref = zi / k0
dscross = []
psi_s = []
for i in range(0, 181):
    psi = 0.0
    y = scs.lpn(ncut,ctheta[i])
    y1 = y[0]
    for n in range(0, ncut +1):
        d_n = (2 * n + 1)
        sum_n = - d_n * y1[n] * tm[n]
        psi = psi + sum_n
    psi = pref * psi
    psi_s = psi_s + [psi]
    dscross = dscross + [psi * scp.conj(psi)]
dscross = scp.real(dscross)


# Generation of the result file:

os.system("del dscross_centered_sphere.txt")
os.system("type nul > dscross_centered_sphere.txt")
fobj1 = open("dscross_centered_sphere.txt", "w")
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
w = scp.imag(psi_s[0])
scat_tot = 4 * scp.pi * w / k0
print()
print("total scattering cross-section: scat_tot = ", scat_tot)


# Plot of "dscross":

if plot_para == 'lg':
    plt.yscale('log')
    plt.plot(theta, dscross, color = 'black', linewidth=2.0)
else:    
    plt.plot(theta, dscross, color = 'black', linewidth=2.0)
plt.xlabel("scattering angle (degree)", fontsize=16)
plt.ylabel("diff. scat. cross-sect.", fontsize=16)
plt.show()
