import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.exp(-0.5*x)*(4 - x) - 2

def df(x):
    return -0.5*(4 - x)*np.exp(-0.5*x) - np.exp(-0.5*x)

plt.clf()
x = np.linspace(0.0, 20.0, 1000)
plt.plot(x, f(x), label="f(x)")
plt.plot(x, df(x), label="df(x)")
plt.grid(True)
plt.savefig("IMG_exe_6_11.pdf")
