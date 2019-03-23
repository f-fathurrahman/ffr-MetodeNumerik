import numpy as np
import matplotlib.pyplot as plt

from chebyshev_poly import *

x = np.array([0.75, 2, 3, 4, 6, 8, 8.5])
y = np.array([1.2, 1.95, 2, 2.4, 2.4, 2.7, 2.6])

coefs = chebyshev_fit(x, y, 4)
print(coefs)

plt.clf()
plt.plot(x, y)
NptsPlot = 200
A = np.min(x)
B = np.max(x)
x_plt = np.linspace(A,B,NptsPlot)
y_plt = np.zeros(NptsPlot)
for i in range(NptsPlot):
    y_plt[i] = chebyshev_fit_eval(coefs, A, B, x_plt[i])
plt.plot(x_plt, y_plt)
plt.savefig("TEMP_test_chebyshev_fit.pdf")


