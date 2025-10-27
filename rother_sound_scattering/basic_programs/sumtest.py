#
# This program provides a test of Eq. (1.80) of the manuscript.
#
import numpy as np
import basics as bas
import math

# Calculation of the separation matrix "S(l,n,nue,x)" according to 
# [Martin, S. 104, Gl. (3.126)]:

def separation_matrix(el, en, nu, x):
    zi = 0 + 1j
    n1 = np.absolute(el) 
    w = zi / 2 / x
    we = np.exp(zi * x) / zi / x
    pf = (-1)**n1 * zi**(en + nu) * we
    j_u = n1
    j_o = en + nu + 1
    S_sumj = 0.0
    for j in range(j_u, j_o):
        du = [0, j - en]
        do = [nu, j - n1]
        s_u = max(du)
        s_o = min(do) + 1
        j1 = math.factorial(j)
        j2 = w**j
        S_sums = 0.        
        for s in range(s_u, s_o):
            jms = j - s
            f1 = np.sqrt(2.0 * nu + 1)
            f1 = np.log(f1)
            f2 = np.sqrt(2. * en + 1)
            f2 = np.log(f2)
            f = f1 + f2
            a1 = 0.
            for i in range(2, nu + s + 1):
                a1 = a1 + np.log(float(i))
            a2 = 0.
            for i in range(2, en + jms + 1):
                a2 = a2 + np.log(float(i))
            a = a1 + a2
            b1 = 0.
            for i in range(2, s + 1):
                b1 = b1 + np.log(float(i))
            b1s = 0.
            for i in range(2, n1 + s + 1):
                b1s = b1s + np.log(float(i))
            b2 = 0.
            for i in range(2, nu - s + 1):
                b2 = b2 + np.log(float(i))
            b3 = 0.
            for i in range(2, jms + 1):
                b3 = b3 + np.log(float(i))
            b3s = 0.
            for i in range(2, jms - n1 + 1):
                b3s = b3s + np.log(float(i))
            b4 = 0.
            for i in range(2, en - jms + 1):
                b4 = b4 + np.log(float(i))
            b = b1s + b1 + b2 + b3s + b3 + b4
            c1 = 0.
            for i in range(2, nu + n1 + 1):
                c1 = c1 + np.log(float(i))                    
            c2 = 0.
            for i in range(2, en - n1 + 1):
                c2 = c2 + np.log(float(i))
            c = c1 + c2
            d1 = 0.
            for i in range(2, nu - n1 + 1):
                d1 = d1 + np.log(float(i))
            d2 = 0.
            for i in range(2, en + n1 + 1):
                d2 = d2 + np.log(float(i))
            d = d1 + d2
            s_fak = np.exp(c - d)
            s_fak = np.sqrt(s_fak)
            s_pref = np.log(s_fak)
            res_fak = np.exp(a + f + s_pref - b)
            S_sums = S_sums + res_fak
        S_sumj = S_sumj + j1 * j2 * S_sums
    S_M = pf * S_sumj
    return S_M    

    
k0b = 1.1 #float(input(' real .................... k0b : '))
n1 = 0 #int(input(' integer .................. n1 : '))
ncut = 10 #int(input(' truncation parameter .... ncut : '))

zi = 0.0 + 1.0j
z_ana = np.exp(zi * k0b) * np.sqrt(4.0 * np.pi * (2.0 * n1 + 1.0))
z_sum = 0.0
for n0 in range(0, ncut + 1):
    a = (-zi)**(n0 + n1)
    b = np.sqrt(4.0 * np.pi * (2.0 * n0 + 1.0))
    #c = bas.SM(0, n0, n1, k0b)
    c = separation_matrix(0, n0, n1, k0b)
    print(f"n0={n0} c = {c}")
    c = np.real(c) #XXX need to be real?
    z_sum += a * b * c

print('analytical value = ', z_ana)
print('approximation = ', z_sum)
print('difference = ', np.abs(z_ana - z_sum))
