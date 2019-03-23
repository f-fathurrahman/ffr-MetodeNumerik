import numpy as np
import matplotlib.pyplot as plt

from chebyshev_poly import *


Npoints = 200
z = np.linspace(-1.0, 1.0, Npoints)


y = np.zeros(Npoints)

plt.clf()
for degree in range(2,7):
    for i in range(Npoints):
        y[i] = chebyshev_poly_eval(degree, z[i])
    plt.plot(z, y, label="degree="+str(degree))

plt.legend()
plt.savefig("TEMP_test_chebyshev_poly_eval.pdf")

