import numpy as np
import matplotlib.pyplot as plt
from newton_interp import *

def func_01(x):
    #return np.sin(2*np.pi*x) - 2*x + 1
    return np.sin(2*np.pi*x)

a = 0.0
b = 1.0
Ninterp = 5
x1 = np.linspace(a, b, Ninterp)
y1 = func_01(x1)

# dense grid
NptsPlot = 100
x2 = np.linspace(a, b, NptsPlot)

coefs = create_newton_polynom(x1, y1)
y2 = np.zeros((NptsPlot,))
for i in range(NptsPlot):
    y2[i] = eval_newton_polynom(coefs, x1, x2[i])

diff = y2 - func_01(x2)

plt.clf()
plt.plot(x1, y1, marker="o", label="sampled")
plt.plot(x2, y2, marker="o", label="interp")
plt.plot(x2, func_01(x2), marker="o", label="true")
plt.grid()
plt.legend()
plt.savefig("TEMP_newton_interp_01.pdf")

plt.clf()
plt.plot(x2, diff, marker="o", label="y2")
plt.grid()
plt.savefig("TEMP_newton_interp_01_diff.pdf")