#            Program "bisphere_z_oriented_input.py"
#

print()
print()
print(" --- Input scattering configuration for acoustic plane wave  ...")
print("   ... scattering on a z-oriented bisphere. Note that the sphere  ...")
print("     ... with radius r=a0 is centered in the laboratoy frame. ---")
print()
print()

import os
os.system("del input_data_bisphere_z_oriented.txt")
os.system("type nul > input_data_bisphere_z_oriented.txt")
fobj = open("input_data_bisphere_z_oriented.txt", "w")


LS = int(input(' 0. iteration (0), 1. iteratio (1), 2. iteration (2) or \
exact solution (3)? '))
print()
para0 = str(input('centered soft (s), hard (h) or penetrable (p) sphere? '))
para1 = str(input('shifted soft (s), hard (h) or penetrable (p) sphere? '))

if (para0 == 'p' and (para1 == 's' or para1 == 'h')):
    k_p_k_0 = float(input('ratio of wave numbers of centered sphere k_p0/k_0: '))
    rho_p_0 = float(input('density of the material inside the centered sphere \
    rho_p0: '))
elif (para1 == 'p' and (para0 == 's' or para0 == 'h')):
    k_p_k_1 = float(input('ratio of wave numbers of shifted sphere k_p1/k_0: '))
    rho_p_1 = float(input('density of the material inside the shifted sphere \
    rho_p1: '))    
elif (para0 == 'p' and para1 == 'p'):
    k_p_k_0 = float(input('ratio of wave numbers of centered sphere k_p0/k_0: '))
    rho_p_0 = float(input('density of the material inside the centered sphere \
    rho_p0: '))
    k_p_k_1 = float(input('ratio of wave numbers of shifted sphere k_p1/k_0: '))
    rho_p_1 = float(input('density of the material inside the shifted sphere \
    rho_p1: '))    
else:
    k_p_k_0 = 1.0
    rho_p_0 = 1.0
    k_p_k_1 = 1.0
    rho_p_1 = 1.0
print()
a0 = float(input('centered sphere with radius a0 [mm]: '))
a1 = float(input('shifted sphere with radius a1 [mm]: '))
beta_0 = float(input('size parameter of centered sphere beta_0: '))
print()
b = float(input('distance between the centers of the spheres b [mm]: '))
print()
ncut = int(input('truncation parameter ... ncut: '))
print()
plot_para = str(input('lin-lin (ll), or lin-log (lg) plot? '))
pl_pm1 = str(input("additional plot of Mie results (y/n)? "))

s00 = str(LS)
s00 = "LS = " + s00 + "\n"
fobj.write(s00)

s0 = "para0 = " + para0 + "\n"
fobj.write(s0)

s1 = "para1 = " + para1 + "\n"
fobj.write(s1)

s2 = str(k_p_k_0)
s2 = "k_p_k_0 = " + s2 + "\n"
fobj.write(s2)

s3 = str(rho_p_0)
s3 = "rho_p_0 = " + s3 + "\n"
fobj.write(s3)

s4 = str(k_p_k_1)
s4 = "k_p_k_1 = " + s4 + "\n"
fobj.write(s4)

s5 = str(rho_p_1)
s5 = "rho_p_1 = " + s5 + "\n"
fobj.write(s5)

s6 = str(a0)
s6 = "a0 = " + s6 + "\n"
fobj.write(s6)

s7 = str(a1)
s7 = "a1 = " + s7 + "\n"
fobj.write(s7)

s8 = str(beta_0)
s8 = "beta_0 = " + s8 + "\n"
fobj.write(s8)

s9 = str(b)
s9 = "b = " + s9 + "\n"
fobj.write(s9)

s10 = str(ncut)
s10 = "ncut = " + s10 + "\n"
fobj.write(s10)

s11 = "plot_para = " + plot_para + "\n"
fobj.write(s11)

s12 = "pl_pm1 = " + pl_pm1
fobj.write(s12)

fobj.close()
