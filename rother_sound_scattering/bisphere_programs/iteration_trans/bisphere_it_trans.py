#            Program "bisphere_it_trans.py"
#
# This program calculates the plane wave scattering behaviour (zero-, first-,  
# and second-order iteration) on a bisphere that is arbitrarily oriented in 
# the laboratory frame. This program is based on a translation of the 
# laboratory frame into the body frame! The configuration is fixed by the 
# distance "b", and the Eulerian angles "\alpha" and "\theta_p". 
#

print()
print()

import numpy as np
import scipy as scp
import scipy.special as scs
import matplotlib.pyplot as plt
import basics as bas
import os
   

# Reading the scattering configuration from input file:

fobj = open("input_data_bisphere_it_trans.txt", "r")

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
alpha = float(z[11])
theta_p = float(z[12])
ncut = int(z[13])
plot_para = z[14]
pl_pm1 = z[15]

k = beta_0 / a0
kb = k * b
beta_1 = k * a1
zi = 0. + 1.0j
pf = - zi / k

theta_l = 360.0
theta_l1 = 361
theta = np.linspace(0.0, theta_l, theta_l1)
ctheta = np.cos(theta * scp.pi / 180.)
ctheta_p = np.cos(theta_p * scp.pi / 180.)
r_theta = theta * scp.pi / 180.


# Calculation of the T-matrices:

if (para0 == 's' and para1 == 's'):
    tm_0 = bas.tm_s(ncut, beta_0)
    tm_1 = bas.tm_s(ncut, beta_1)
elif (para0 == 's' and para1 == 'h'):
    tm_0 = bas.tm_s(ncut, beta_0)
    tm_1 = bas.tm_h(ncut, beta_1)
elif (para0 == 's' and para1 == 'p'):
    tm_0 = bas.tm_s(ncut, beta_0)
    tm_1 = bas.tm_p(ncut, beta_1, a1, k, k_p_k_1, rho_p_1)
elif (para0 == 'h' and para1 == 's'):
    tm_0 = bas.tm_h(ncut, beta_0)
    tm_1 = bas.tm_s(ncut, beta_1)
elif (para0 == 'h' and para1 == 'h'):
    tm_0 = bas.tm_h(ncut, beta_0)
    tm_1 = bas.tm_h(ncut, beta_1)
elif (para0 == 'h' and para1 == 'p'):
    tm_0 = bas.tm_h(ncut, beta_0)
    tm_1 = bas.tm_p(ncut, beta_1, a1, k, k_p_k_1, rho_p_1)
elif (para0 == 'p' and para1 == 's'):
    tm_0 = bas.tm_p(ncut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_s(ncut, beta_1)
elif (para0 == 'p' and para1 == 'h'):
    tm_0 = bas.tm_p(ncut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_h(ncut, beta_1)
else:
    tm_0 = bas.tm_p(ncut, beta_0, a0, k, k_p_k_0, rho_p_0)
    tm_1 = bas.tm_p(ncut, beta_1, a1, k, k_p_k_1, rho_p_1)
t_o = tm_0
t_o1 = tm_1


# Calculation of the separation matrix:

sm_l_i_j = []   
for l in range(0, ncut + 1):   
    sm_i_j = []
    for i in range(0, ncut + 1):
        sm_j = []
        for j in range(0, ncut + 1):
            sm_j = sm_j + [bas.SM(l, i, j, kb)]
        sm_i_j = sm_i_j + [sm_j]
    sm_l_i_j = sm_l_i_j + [sm_i_j]


# Mie coefficients of the centered and shifted sphere:
    
c_m0_mie = []
q_m0_mie = []   
for m0 in range(0, ncut + 1):
    d_m0 = np.sqrt(4. * np.pi * (2 * m0 + 1)) * zi**m0
    c_m0_mie = c_m0_mie + [t_o[m0] * d_m0]
    q_m0_mie = q_m0_mie + [np.exp(zi * kb * ctheta_p) * t_o1[m0] *\
                           d_m0]


# Calculation of the differential scattering cross-section "dscross" 
# as a function of the scattering angle "theta" in the interval [0, 2 \pi]
# in steps of 1 degree: only Mie theory for intercomparison purposes!

dscross_0 = []
dscross_1 = []
for i in range(0, theta_l1):
    psi_0_0 = 0.0
    psi_1_0 = 0.0
    
    for m0 in range(0, ncut + 1):
        Y_0m0 = scs.sph_harm(0,m0,0.,r_theta[i])
        psi_0_0 = psi_0_0 + c_m0_mie[m0] * (-zi)**m0 * Y_0m0
    psi_0_0 = pf * psi_0_0
    dscross_0 = dscross_0 + [psi_0_0 * np.conj(psi_0_0)] 
                             
    for m1 in range(0, ncut + 1):
        Y_0m1 = scs.sph_harm(0,m1,0.,r_theta[i])
        psi_1_0 = psi_1_0 + q_m0_mie[m1] * (-zi)**m1 * Y_0m1
    psi_1_0 = pf * psi_1_0
    dscross_1 = dscross_1 + [psi_1_0 * np.conj(psi_1_0)]
    
    
# Iterative solutions in each local system:
    
c_l0_m0_0 = []
q_l1_n1_0 = []
for l0 in range(-ncut, ncut + 1):    
    c_m0_1 = []
    q_n1_1 = []
    for m0 in range(np.abs(l0), ncut + 1):
        c_m0_0 = (c_m0_mie[m0] if l0 == 0 else 0.0)
        q_n1_0 = (q_m0_mie[m0] if l0 == 0 else 0.0)
        c_m0_1 = c_m0_1 + [c_m0_0]
        q_n1_1 = q_n1_1 + [q_n1_0]
    c_l0_m0_0 = c_l0_m0_0 + [c_m0_1]
    q_l1_n1_0 = q_l1_n1_0 + [q_n1_1]

c_l0_m0_1 = bas.Q12_tl_it1(ncut, alpha, theta_p, sm_l_i_j, t_o, q_l1_n1_0)    
q_l1_n1_1 = bas.C21_tl_it1(ncut, alpha, theta_p, sm_l_i_j, t_o1, c_l0_m0_0) 


c_it_1 = []
for l in range(-ncut, ncut + 1):
    v1 = c_l0_m0_0[l + ncut]
    v2 = c_l0_m0_1[l + ncut]
    v_it_1 = np.array(v1) - np.array(v2)
    c_it_1 = c_it_1 + [v_it_1]
q_it_1 = []
for l in range(-ncut, ncut + 1):
    w1 = q_l1_n1_0[l + ncut]
    w2 = q_l1_n1_1[l + ncut]
    w_it_1 = np.array(w1) - np.array(w2)
    q_it_1 = q_it_1 + [w_it_1]
       
cc1 = bas.C21_tl_it1(ncut, alpha, theta_p, sm_l_i_j, t_o1, c_it_1)
qq1 = bas.Q12_tl_it1(ncut, alpha, theta_p, sm_l_i_j, t_o, q_it_1)

c_l0_m0_2 = bas.Q12_tl_it1(ncut, alpha, theta_p, sm_l_i_j, t_o, cc1)    
q_l1_n1_2 = bas.C21_tl_it1(ncut, alpha, theta_p, sm_l_i_j, t_o1, qq1)

c_final = []
for l0 in range(-ncut, ncut + 1):
    v1 = c_l0_m0_0[l0 + ncut]
    v2 = c_l0_m0_1[l0 + ncut]
    v3 = c_l0_m0_2[l0 + ncut]
    if LS == 0:
        v = v1
    elif LS == 1:
        v = np.array(v1) - np.array(v2)
    else:
        v = np.array(v1) - np.array(v2) + np.array(v3)
    c_final = c_final + [v]

q_final = []
for l1 in range(-ncut, ncut + 1):
    w1 = q_l1_n1_0[l1 + ncut]
    w2 = q_l1_n1_1[l1 + ncut]
    w3 = q_l1_n1_2[l1 + ncut]
    if LS == 0:
        w = w1
    elif LS == 1:
        w = np.array(w1) - np.array(w2)
    else:
        w = np.array(w1) - np.array(w2) + np.array(w3)
    q_final = q_final + [w]


# Calculation of the differential scattering cross-section "dscross_g" 
# as a function of the scattering angle "theta" in the interval [0, 2 \pi]
# in steps of 1.0 degree (in the laboratory frame):
   
psi_s_g = []
dscross_g = []
for i in range(0, theta_l1):
    c_k = kb * (np.cos(r_theta[i]) * ctheta_p + np.sin(r_theta[i]) * \
    np.sin(theta_p * scp.pi / 180.) * np.cos(alpha *scp.pi / 180.)) 
    psi_s = 0.0    
    for l0 in range(-ncut, ncut + 1):
        c_l0 = c_final[l0 + ncut]
        q_l0 = q_final[l0 + ncut]
        for m0 in range(np.abs(l0), ncut + 1):
            m0z = m0 - np.abs(l0)
            pref = (-zi)**m0
            if i <= 180:
                Y_l0_m0 = scs.sph_harm(l0,m0,0.,r_theta[i])
            else:
                Y_l0_m0 = scs.sph_harm(l0,m0,np.pi,r_theta[360 - i])
            psi_s = psi_s + (c_l0[m0z] + np.exp(-zi * c_k) * q_l0[m0z]) * \
                             pref * Y_l0_m0
    psi_s = pf * psi_s
    psi_s_g = psi_s_g + [psi_s]
    dscross_g = dscross_g + [psi_s * np.conj(psi_s)]                             
dscross_g = np.real(dscross_g)   
dscross_0 = np.real(dscross_0)
dscross_1 = np.real(dscross_1)
dscross_zw = 2. * (dscross_0 + dscross_1)
dscross_zw1 = dscross_0 + dscross_1   #  für a1 << a0 !
        

# Generation of the result file:

os.system("del dscross_bisphere_it_trans.txt")
os.system("type nul > dscross_bisphere_it_trans.txt")
fobj1 = open("dscross_bisphere_it_trans.txt", "w")
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
        plt.plot(theta, dscross_zw, '--', color = 'black', linewidth=2.0, \
                 label = "2 * (Mie_a0 + Mie_a1)")
        plt.plot(theta, dscross_g, color = 'black', linewidth=2.0, \
                 label = "bisphere")
    else:
        plt.plot(theta, dscross_g, color = 'black', linewidth=2.0, \
                 label = "bisphere")
else:
    if pl_pm1 == 'y':
        plt.plot(theta, dscross_zw, '--', color = 'black', linewidth=2.0, \
                 label = "2 * (Mie_a0 + Mie_a1)")
        plt.plot(theta, dscross_g, color = 'black', linewidth=2.0, \
                 label = "bisphere")
    else:
        plt.plot(theta, dscross_g, color = 'black', linewidth=2.0, \
                 label = "bisphere")
plt.xlabel("scattering angle [deg.]", fontsize=16)
plt.ylabel("diff. scat. cross-sect.", fontsize=16)
plt.legend()
#plt.legend(loc = 'lower right')
plt.legend(loc = 'upper center')
#plt.legend(loc = 'upper left')
plt.show()
 