import numpy as np
from polyN_fit import *

x = np.array([0.0, 1.0, 2.0, 2.5, 3.0])
y = np.array([2.9, 3.7, 4.1, 4.4, 5.0])

coefs = polyN_fit(x, y, 1)
print(coefs)

import matplotlib.pyplot as plt

plt.clf()
plt.plot(x, y, linewidth=0, marker="o", label="data")
xplt = np.linspace(x[0], x[-1], 500)
plt.plot(xplt, np.polyval(np.poly1d(coefs[::-1]), xplt), label="linreg" )
plt.legend()
plt.savefig("TEMP_test_polyN_fit_01.pdf")
