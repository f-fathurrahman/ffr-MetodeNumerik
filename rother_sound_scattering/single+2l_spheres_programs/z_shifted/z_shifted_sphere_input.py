#            Program "z_shifted_sphere_input.py"
#

print()
print('        --- Input scattering configuration for acoustic plane ...')
print('            ... wave scattering on a z-shifted sphere! ---')
print()

import os
os.system("del input_data_z_shifted_sphere.txt")
os.system("type nul > input_data_z_shifted_sphere.txt")
fobj = open("input_data_z_shifted_sphere.txt", "w")

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
b = float(input('shift along the z-axis b>=2a [mm]: '))
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

s5 = str(b)
s5 = "b = " + s5 + "\n"
fobj.write(s5)

s6 = str(ncut)
s6 = "ncut = " + s6 + "\n"
fobj.write(s6)

s7 = "plot_para = " + plot_para
fobj.write(s7)

fobj.close()
