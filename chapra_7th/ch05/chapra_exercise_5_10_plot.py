import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**4 - 8*x**3 - 35*x**2 + 450*x - 1001

x = np.linspace(0.0, 10.0, 1000)
plt.clf()
plt.plot(x, f(x))
plt.grid(True)
plt.savefig("IMG_exe_5_10.pdf")

x = np.linspace(4.0, 6.0, 1000)
plt.clf()
plt.plot(x, f(x))
plt.grid(True)
plt.savefig("IMG_exe_5_10_v2.pdf")