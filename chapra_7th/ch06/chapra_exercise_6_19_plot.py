import numpy as np
import matplotlib.pyplot as plt

R = 3.0  # m
V = 30.0 # m^3

def f(h):
    return V - np.pi*h**2*(3*R - h)/3

h = np.linspace(0.0, 10.0, 1000)
plt.clf()
plt.plot(h, f(h))
plt.grid(True)
plt.savefig("IMG_exe_6_19.pdf")