import numpy as np
import matplotlib.pyplot as plt

r = 1.0
ρs = 200.0
ρw = 1000.0

def f(h):
    Vw = 4*np.pi*r**3/3 - np.pi*h**2/3*(3*r - h) # displaced volume of water
    Vs = 4*np.pi*r**3/3
    return ρw*Vw - ρs*Vs

h = np.linspace(0.0, 2*r, 1000)
plt.clf()
plt.plot(h, f(h))
plt.grid(True)
plt.savefig("IMG_exe_5_19.pdf")