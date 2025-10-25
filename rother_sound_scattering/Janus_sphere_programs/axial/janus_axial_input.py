#            Program "janus_axial_input.py"
print()
print()
print(" --- Acoustic plane wave scattering on a Janus sphere! ---")
print("           ... (axisymmetric orientation) ---")
print()
print()
print()

import os
os.system("del input_data_janus_axial.txt")
os.system("type nul > input_data_janus_axial.txt")
fobj = open("input_data_janus_axial.txt", "w")


# Input scattering configuration:

para = str(input('hard-penetrable (hp), soft-penetrable (sp), or \
hard-soft (hs) Janus sphere? '))
print()
if para == 'hp' or para == 'sp':
    k_p_k = float(input('ratio of wave numbers k_p/k_0: '))
    rho_p = float(input('density of the material inside the sphere rho_p: '))
else:
    k_p_k = 1.0
    rho_p = 1.0
a = float(input('radius ... a [mm]: '))
beta = float(input('size parameter ... beta: '))
theta_j = float(input('splitting angle ... theta_j [degree]: '))
print()
ncut = int(input("truncation parameter ...  ncut: "))
print()
plot_para = str(input('lin-lin (ll), or lin-log (lg) plot? '))

s0 = "para = " + para + "\n"
fobj.write(s0)

s1 = str(k_p_k)
s1 = "k_p_k = " + s1 + "\n"
fobj.write(s1)

s2 = str(rho_p)
s2 = "rho_p = " + s2 + "\n"
fobj.write(s2)

s3 = str(a)
s3 = "a = " + s3 + "\n"
fobj.write(s3)

s4 = str(beta)
s4 = "beta = " + s4 + "\n"
fobj.write(s4)

s5 = str(theta_j)
s5 = "theta_j = " + s5 + "\n"
fobj.write(s5)

s6 = str(ncut)
s6 = "ncut = " + s6 + "\n"
fobj.write(s6)

s7 = "plot_para = " + plot_para
fobj.write(s7)

fobj.close()
