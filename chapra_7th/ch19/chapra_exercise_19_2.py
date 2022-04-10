import numpy as np
import matplotlib.pyplot as plt

from fourier_fit import *

x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=np.float64)
y = np.array([114, 188, 245, 311, 351, 359, 308, 287, 260, 211, 159, 131], dtype=np.float64)

print("Ndata = ", len(x))

A0, A1, B1, omega = sinusoid_fit(x, y)

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
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
plt.xticks(np.arange(1,13), months)
plt.legend()
plt.grid()
plt.savefig("IMG_chapra_exercise_19_2.pdf")


print("Pada pertengahan Agustus:")
print(sinusoid_fit_eval(A0, A1, B1, omega, 8.5))

