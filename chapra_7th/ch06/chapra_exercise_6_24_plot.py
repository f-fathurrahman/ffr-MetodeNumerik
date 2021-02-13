import numpy as np
import matplotlib.pyplot as plt

def f1(x):
    return x**2 + 1

def f2(x):
    return 2*np.cos(x)

x = np.linspace(-2.0, 2.0, 1000)
y1 = f1(x)
y2 = f2(x)
plt.plot(x, y1, label="f1")
plt.plot(x, y2, label="f2")
plt.grid(True)
plt.savefig("IMG_exe_6_24.pdf")