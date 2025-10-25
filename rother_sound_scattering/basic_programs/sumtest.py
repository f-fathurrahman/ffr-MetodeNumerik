#            Program "sumtest.py"
#
# This program provides a test of Eq. (1.80) of the manuscript.
#

print()
print()
print("                      --- Sum test! ---")
print()
print()

import numpy as np
import basics as bas

    
k0b = 1.1 #float(input(' real .................... k0b : '))
n1 = 3 #int(input(' integer .................. n1 : '))
ncut = 10 #int(input(' truncation parameter .... ncut : '))

zi = 0. + 1.0j
z_ana = np.exp(zi * k0b) * np.sqrt(4. * np.pi * (2. * n1 + 1.))
z_sum = 0.
for n0 in range(0, ncut + 1):
    a = (-zi)**(n0 + n1)
    b = np.sqrt(4. * np.pi * (2. * n0 + 1.))
    c = bas.SM(0, n0, n1, k0b)
    c = np.real(c)
    z_sum = z_sum + a * b * c
print('analytical value:   ', z_ana)
print('approximation:      ', z_sum)
