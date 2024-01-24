import numpy as np
import matplotlib.pyplot as plt

from fourier_fit import *

x = np.array([0.5, 0.7, 0.9, 1.1, 1.3])
y = np.array([6.8, 3.2, -4.1, -3.9, 3.3])

A0, A, B = fourier_fit(x, y, m=1)

print("A0 = ", A0)
print("A  = ", A)
print("B  = ", B)


xmin = np.min(x)
xmax = np.max(x)
NptsPlt = 200
x_plt = np.linspace(xmin, xmax, NptsPlt)
y_plt = np.zeros(NptsPlt)
for i in range(NptsPlt):
    y_plt[i] = fourier_fit_eval(A0, A, B, x, x_plt[i])

plt.clf()
plt.plot(x, y, "o", label="data")
plt.plot(x_plt, y_plt, label="fourier-fit")
plt.legend()
plt.savefig("TEMP_test_fourier_fit.pdf")

