#            Program "plot_Bessel_function"
#
# This program plots the spherical Bessel function for real-valued arguments 
# z in the interval [0.,15.]. The results can be compared to Fig. 10.1 
# in reference [6].
#

print()
print()
print("           --- Plot of spherical Bessel functions! ---")
print()
print()

import matplotlib.pyplot as plt
import numpy as np
import scipy.special as scs


# Input order of Bessel function:

n = int(input('Order of Bessel function ... n: '))
zi = 0. + 1.0j

# Calculation of the Bessel function of order "n" in [0,14]:
    
z = np.linspace(0.0, 150.0, 1501)
bes = scs.spherical_jn(n,z,derivative=False)

# Plot of Bessel function:
    
#plt.yscale('log')  # use this command for lin-log plot!
plt.plot(z, bes, color = 'red', linewidth=2.0)
plt.xlabel("z", fontsize=16)
plt.ylabel("spherical Bessel function j_n(z)", fontsize=16)
plt.show()

