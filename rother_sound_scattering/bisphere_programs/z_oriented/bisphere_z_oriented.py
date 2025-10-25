#          Program "bisphere_z_oriented.py"
#
# Acoustic plane wave scattering on a z-oriented bisphere. Note that the 
# sphere with radius "r=a0" is centered in the laboratory frame. 
#

print()
print()

import numpy as np
import scipy as scp
import scipy.special as scs
import scipy.linalg as scl
import matplotlib.pyplot as plt
import basics as bas
import os


# Reading scattering configuration from input file:

fobj = open("input_data_bisphere_z_oriented.txt", "r")

z = []
for line in fobj:
    line = line.strip()
    arr = line.split("= ")
    wert = str(arr[1])
    z = z + [wert]
fobj.close()

LS = int(z[0])
para0 = z[1]
para1 = z[2]
k_p_k_0 = float(z[3])
rho_p_0 = float(z[4])
k_p_k_1 = float(z[5])
rho_p_1 = float(z[6])
a0 = float(z[7])
a1 = float(z[8])
beta_0 = float(z[9])
b = float(z[10])
ncut = int(z[11])
plot_para = z[12]
pl_pm1 = z[13]

n1cut = ncut
theta_p = 0.0
zi = 0. + 1.0j
k = beta_0 / a0
kb = k * b
beta_1 = k * a1


# Calculation of the T-matrices:

if (para0 == 's' and para1 == 's'):
    tm_0 = bas.tm_s(ncut, beta_0)
    tm_1 = bas.tm_s(n1cut, beta_1)
elif (para0 == 's' and para1 == 'h'):
    tm_0 = bas.tm_s(ncut, beta_0)
    tm_1 = bas.tm_h(n1cut, beta_1)
elif (para0 == 's' and para1 == 'p'):
    tm_0 = bas.tm_s(ncut, beta_0)
    tm_1 = bas.tm_p(n1cut, beta_1, a1, k, k_p_k_1, rho_p_1)
elif (para0 == 'h' and para1 == 's'):
    tm_0 = bas.tm_h(ncut, beta_0)
    tm_1 = bas.tm_s(n1cut, beta_1)
elif (para0 == 'h' and para1 == 'h'):
    tm_0 = bas.tm_h(ncut, beta_0)
    tm_1 = bas.tm_h(n1cut, beta_1)
elif (para0 == 'h' and para1 == 'p'):
    tm_0 = bas.tm_h(ncut, beta_0)
    tm_1 = bas.tm_p(n1cut, beta_1, a1, k, k_p_k_1, rho_p_1)
elif (para0 == 'p' and para1 == 's'):
    tm_0 = bas.tm_p(ncut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_s(n1cut, beta_1)
elif (para0 == 'p' and para1 == 'h'):
    tm_0 = bas.tm_p(ncut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_h(n1cut, beta_1)
else:
    tm_0 = bas.tm_p(ncut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_p(n1cut, beta_1, a1, k, k_p_k_1, rho_p_1)


# Calculation of the separation matrix:
    
y_m0_m1 = []   
for i in range(0, ncut + 1):   
    y_os = []
    for j in range(0, ncut + 1):
        y_os = y_os + [bas.SM(0, i, j, kb)]
    y_m0_m1 = y_m0_m1 + [y_os]


# Mie calculation for single spheres:
    
c_m0_0 = []  
for m0 in range(0, ncut + 1):
    d_m0 = np.sqrt(4. * np.pi * (2 * m0 + 1)) * zi**m0
    c = tm_0[m0] * d_m0
    c_m0_0 = c_m0_0 + [c]

q_m1_0 = []
for m1 in range(0, n1cut + 1):
    d_m1 = np.sqrt(4. * np.pi * (2 * m1 + 1)) * zi**m1
    q = tm_1[m1] * d_m1
    q_m1_0 = q_m1_0 + [q]


# Calculation of the matrices of the T-matrix equation:
    
c11 = []
for ic11 in range(0, ncut + 1):
    c1 = []
    for jc11 in range(0, ncut + 1):
        c = (1.0 if jc11 == ic11 else 0.0)
        c1 = c1 + [c]
    c11 = c11 + [c1]
    
q22 = []
for iq22 in range(0, n1cut + 1):
    q2 = []
    for jq22 in range(0, n1cut + 1):
        q = (1.0 if jq22 == iq22 else 0.0)
        q2 = q2 + [q]
    q22 = q22 + [q2]
    
q12 = []
c0 = []
for iq12 in range(0, ncut + 1):
    q1 = []
    S_m0 = y_m0_m1[iq12]
    q_to = - tm_0[iq12] * (-1)**iq12
    d_ic = np.sqrt(4. * np.pi * (2 * iq12 + 1)) * zi**iq12
    c0 = c0 + [tm_0[iq12] * d_ic]
    for jq12 in range(0, n1cut + 1):
        q = q_to * S_m0[jq12]
        q1 = q1 + [q]
    q12 = q12 + [q1]
     
c21 = []
q0 = []
for ic21 in range(0, n1cut + 1):
    c2 = []
    S_m0 = y_m0_m1[ic21]
    c_to1 = - tm_1[ic21]
    d_iq = np.sqrt(4. * np.pi * (2 * ic21 + 1)) * zi**ic21
    q0 = q0 + [np.exp(zi * kb) * tm_1[ic21] * d_iq]
    for jc21 in range(0, ncut + 1):
        c = c_to1 * S_m0[jc21] * (-1)**jc21
        c2 = c2 + [c]
    c21 = c21 + [c2]

c11 = np.array(c11)
q22 = np.array(q22)
c21 = np.array(c21)
q12 = np.array(q12)
c0 = np.array(c0)
q0 = np.array(q0)

cq = np.matmul(c21, q12)
qc = np.matmul(q12, c21)


# Solving the T-matrix equation (iteration and rigorous):
    
if LS == 0:
    mc = c11
    mq = q22
    ci = c0
    qi = q0
    c_m0_00 = scl.solve(mc,ci)
    q_m1_00 = scl.solve(mq,qi)
elif LS == 1:
    mc = c11
    mq = q22
    ci = c0 - q12.dot(q0)
    qi = q0 - c21.dot(c0)    
    c_m0_00 = scl.solve(mc,ci)
    q_m1_00 = scl.solve(mq,qi)
elif LS == 2:
    mc = c11
    mq = q22
    ci1 = c0 - q12.dot(q0)
    qi1 = q0 - c21.dot(c0)    
    ci_zw = np.matmul(c21, ci1)
    ci = ci1 + np.matmul(q12, ci_zw)
    qi_zw = np.matmul(q12, qi1)
    qi = qi1 + np.matmul(c21, qi_zw)
    c_m0_00 = scl.solve(mc,ci)
    q_m1_00 = scl.solve(mq,qi)    
elif LS == 3:
    mc = c11 - qc
    mq = q22 - cq
    ci = c0 - q12.dot(q0)
    qi = q0 - c21.dot(c0)
    c_m0_00 = scl.solve(mc,ci)
    q_m1_00 = scl.solve(mq,qi)        
c_m0 = c_m0_00
q_m1 = q_m1_00
               

# Calculation of the differential scattering cross-section "dscross_g" 
# as a function of the scattering angle "theta" in the interval [0, \pi]
# in steps of 1 degree. 

theta_l = 180.0
theta_l1 = 181
theta = np.linspace(0.0, theta_l, theta_l1)
r_theta = theta * scp.pi / 180.
ctheta = np.cos(r_theta)
zi = 0. + 1.0j
pref = - zi / k
psi_s_g = []
dscross_g = []
dscross_0 = []
dscross_1 = []
for i in range(0, 181):
    psi_0 = 0.0
    psi_0_1 = 0.0
    psi_1 = 0.0
    psi_0_0 = 0.0
    
    for m0 in range(0, ncut + 1):
        Y_0m0 = scs.sph_harm(0,m0,0.,r_theta[i])
        psi_0 = psi_0 + c_m0[m0] * (-zi)**m0 * Y_0m0
        psi_0_0 = psi_0_0 + c_m0_0[m0] * (-zi)**m0 * Y_0m0
    psi_0 = pref * psi_0
    psi_0_0 = pref * psi_0_0
    dscross_0 = dscross_0 + [psi_0_0 * np.conj(psi_0_0)] 
                             
    for m1 in range(0, n1cut + 1):
        Y_0m1 = scs.sph_harm(0,m1,0.,r_theta[i])
        psi_0_1 = psi_0_1 + q_m1_0[m1] * (-zi)**m1 * Y_0m1
    psi_0_1 = pref * psi_0_1
    dscross_1 = dscross_1 + [psi_0_1 * np.conj(psi_0_1)]

    for m1s in range(0, n1cut + 1):
        Y_0m1s = scs.sph_harm(0,m1s,0.,r_theta[i])
        sum_m1s = q_m1[m1s] * (-zi)**m1s * Y_0m1s 
        psi_1 = psi_1 + sum_m1s
    c_k = kb * np.cos(np.abs(theta_p - theta[i]) * scp.pi / 180.)
    psi_1 = pref * psi_1 * np.exp(- zi * c_k)
    psi_s = psi_0 + psi_1
    psi_s_g = psi_s_g + [psi_s]
    dscross_g = dscross_g + [psi_s * np.conj(psi_s)]

dscross_g = np.real(dscross_g)
dscross_0 = np.real(dscross_0)
dscross_1 = np.real(dscross_1)
dscross_zw = 2. * (dscross_0 + dscross_1)
dscross_zw1 = dscross_0 + dscross_1   #  für a1 << a0 !


# Generation of the result file:

os.system("del dscross_bisphere_z_oriented.txt")
os.system("type nul > dscross_bisphere_z_oriented.txt")
fobj1 = open("dscross_bisphere_z_oriented.txt", "w")
nl = "\n"
for i in range(0, theta_l1):
    a1 = str(theta[i])
    a2 = " = "
    a3 = str(dscross_g[i]) + nl
    a = a1 + a2 + a3
    fobj1.write(a)
fobj1.close()


# Calculation of the total scattering cross-section "scat_tot" by use 
# of the optical theorem:

print()
print()
print('Results: ')
w = np.imag(psi_s_g[0])
scat_tot = 4 * scp.pi * w / k
print()
print("total scattering cross-section: scat_tot = ", scat_tot)
print()


# Plot of the results:
    
if plot_para == 'lg':
    plt.yscale('log')
    if pl_pm1 == 'y':
        plt.plot(theta, dscross_zw, '-.', color = 'black', linewidth=2.0, \
                 label = "2 * (Mie_a0 + Mie_a1)")
        plt.plot(theta, dscross_g, color='black', linewidth=2.0, \
                 label = "bisphere")
    else:
        plt.plot(theta, dscross_g, linewidth=2.0, label = "bisphere")
else:
    if pl_pm1 == 'y':
        plt.plot(theta, dscross_zw, '-.', color = 'black', linewidth=2.0, \
                 label = "2 * (Mie_a0 + Mie_a1)")
        plt.plot(theta, dscross_g, color='black', linewidth=2.0, \
                 label = "bisphere")
    else:
        plt.plot(theta, dscross_g, color = 'black', linewidth=2.0, \
                 label = "bisphere")
plt.xlabel("scattering angle [degree]", fontsize=16)
plt.ylabel("diff. scat. cross-sect.", fontsize=16)
plt.legend(loc = 'upper center')
# plt.legend()
plt.show()

