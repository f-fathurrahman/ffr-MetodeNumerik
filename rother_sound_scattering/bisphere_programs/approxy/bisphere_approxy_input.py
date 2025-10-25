#            Program "bisphere_approxy_input.py"
#

print()
print()
print('     --- Input scattering configuration for acoustic plane ...')
print('      ... wave scattering on combinations of acoustically ...')
print("      ... soft, hard, or penetrable spheres arbitrarily ...")
print("      ... oriented in the laboratory frame (no interaction is ...")
print(" ...                       considered)! ---")
print()
print()

import os
os.system("del input_data_bisphere_approxy.txt")
os.system("type nul > input_data_bisphere_approxy.txt")
fobj = open("input_data_bisphere_approxy.txt", "w")

para0 = str(input('centered soft (s), hard (h) or penetrable (p) sphere? '))
para1 = str(input('shifted soft (s), hard (h) or penetrable (p) sphere? '))

if (para0 == 'p' and (para1 == 's' or para1 == 'h')):
    k_p_k_0 = float(input('ratio of wave numbers of centered sphere k_p0/k: '))
    rho_p_0 = float(input('density of the material inside the centered sphere \
    rho_p0: '))
elif (para1 == 'p' and (para0 == 's' or para0 == 'h')):
    k_p_k_1 = float(input('ratio of wave numbers of shifted sphere k_p1/k: '))
    rho_p_1 = float(input('density of the material inside the shifted sphere \
    rho_p1: '))    
elif (para0 == 'p' and para1 == 'p'):
    k_p_k_0 = float(input('ratio of wave numbers of centered sphere k_p0/k: '))
    rho_p_0 = float(input('density of the material inside the centered sphere \
    rho_p0: '))
    k_p_k_1 = float(input('ratio of wave numbers of shifted sphere k_p1/k: '))
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
alpha = float(input('Eulerian angle of orientation alpha [deg.]: '))
theta_p = float(input('Eulerian angle of orientation theta_p [deg]: '))
print()
ppm = int(input('dif. scat. cross. between [0,pi] (0) or [0,2*pi] (1)? '))
plot_para = str(input('lin-lin (ll), or lin-log (lg) plot? '))

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

s10 = str(alpha)
s10 = "alpha = " + s10 + "\n"
fobj.write(s10)

s11 = str(theta_p)
s11 = "theta_p = " + s11 + "\n"
fobj.write(s11)

s12 = str(ppm)
s12= "ppm = " + s12 + "\n"
fobj.write(s12)

s13 = "plot_para = " + plot_para
fobj.write(s13)

fobj.close()
