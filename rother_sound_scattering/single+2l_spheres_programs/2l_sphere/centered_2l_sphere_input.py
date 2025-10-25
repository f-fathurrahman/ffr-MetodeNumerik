#            Program "centered_2l_sphere_input.py"
#

print()
print()
print('    --- Input scattering configuration for acoustic plane ...')
print('       ... wave scattering on a centered 2layer sphere! ---')
print()
print()
print()

import os
os.system("del input_data_centered_2l_sphere.txt")
os.system("type nul > input_data_centered_2l_sphere.txt")
fobj = open("input_data_centered_2l_sphere.txt", "w")

# Input scattering configuration:

para = str(input('soft+penetrable (sp) or hard+penetrable (hp)? '))
print()
k_p_k = float(input('ratio of wave numbers k_p/k0: '))
rho_p = float(input('density of the material inside the sphere ... rho_p: '))
a = float(input('radius of the core ... a [mm]: '))
ts = float(input('thickness of the shell ... ts [mm]: '))
beta_b = float(input('size parameter of the 2l-sphere ... beta_b: '))
print()
ncut = int(input('truncation parameter  ... ncut: '))
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

s4 = str(ts)
s4 = "ts = " + s4 + "\n"
fobj.write(s4)

s5 = str(beta_b)
s5 = "beta_b = " + s5 + "\n"
fobj.write(s5)

s6 = str(ncut)
s6 = "ncut = " + s6 + "\n"
fobj.write(s6)

s7 = "plot_para = " + plot_para
fobj.write(s7)

fobj.close()
