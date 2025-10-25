#            Program "centered_sphere_input.py"
#

print()
print()
print('     --- Input scattering configuration for acoustic plane ...')
print('            ... wave scattering on a centered sphere! ---')
print()
print()

import os
os.system("del input_data_centered_sphere.txt")
os.system("type nul > input_data_centered_sphere.txt")
fobj = open("input_data_centered_sphere.txt", "w")

para = str(input('soft (s), hard (h), or penetrable (p) sphere? '))
print()
if para == 'p':
    k_p_k = float(input('ratio of wave numbers k_p/k_0: '))
    rho_p = float(input('density of the material inside the sphere rho_p: '))
else:
    k_p_k = 1.0
    rho_p = 1.0
a = float(input('radius a [mm]: '))
beta = float(input('size parameter beta: '))
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

s5 = "plot_para = " + plot_para
fobj.write(s5)

fobj.close()
