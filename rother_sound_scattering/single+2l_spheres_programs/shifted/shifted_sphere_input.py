#            Program "shifted_sphere_input.py"
#

print()
print('     --- Input scattering configuration for acoustic plane ...')
print('            ... wave scattering on a shifted sphere! ---')
print()

import os
os.system("del input_data_shifted_sphere.txt")
os.system("type nul > input_data_shifted_sphere.txt")
fobj = open("input_data_shifted_sphere.txt", "w")

para = str(input('soft (s), hard (h), or penetrable (p) sphere? '))
print()
if para == 'p':
    k_p_k = float(input('ratio of wave numbers k_p/k0: '))
    rho_p = float(input('density of the material inside the sphere rho_p: '))
else:
    k_p_k = 1.0
    rho_p = 1.0
a = float(input('radius a [mm]: '))
beta = float(input('size parameter beta: '))
print()
alpha = float(input('Eulerian angle alpha [degree]: '))
theta_p = float(input('Eulerian angle theta_p [degree]: '))
b = float(input('shift along new z-axis  b>0 [mm]: '))
print()
ncut = int(input('truncation parameter ncut: '))
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

s5 = str(alpha)
s5 = "alpha = " + s5 + "\n"
fobj.write(s5)

s6 = str(theta_p)
s6 = "theta_p = " + s6 + "\n"
fobj.write(s6)

s7 = str(b)
s7 = "b = " + s7 + "\n"
fobj.write(s7)

s8 = str(ncut)
s8 = "ncut = " + s8 + "\n"
fobj.write(s8)

s9 = "plot_para = " + plot_para
fobj.write(s9)

fobj.close()
