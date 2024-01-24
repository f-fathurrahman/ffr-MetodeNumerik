import numpy as np
import matplotlib.pyplot as plt

from fourier_fit import *

x = np.array([0, 2, 4, 5, 7, 9, 12, 15, 20, 22, 24], dtype=np.float64)
y = np.array([7.6, 7, 7.1, 6.5, 7.4, 7.2, 8.9, 8.8, 8.9, 7.9, 7.0], dtype=np.float64)

print("Ndata = ", len(x))

A0, A1, B1, omega = sinusoid_fit(x, y, equispaced=False, T=24.0)

print("A0 = ", A0)
print("A1 = ", A1)
print("B1 = ", B1)
print("omega = ", omega)

xmin = np.min(x)
xmax = np.max(x)
NptsPlt = 200
x_plt = np.linspace(xmin, xmax, NptsPlt)
y_plt = np.zeros(NptsPlt)
for i in range(NptsPlt):
    y_plt[i] = sinusoid_fit_eval(A0, A1, B1, omega, x_plt[i])

plt.clf()
plt.plot(x, y, "o", label="data")
plt.plot(x_plt, y_plt, label="sinusoid-fit")
plt.xticks(np.arange(0,25,2))
plt.grid()
plt.legend()
plt.savefig("TEMP_chapra_19_3.pdf")

