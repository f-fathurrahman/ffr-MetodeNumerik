import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.sin(10*x) + np.cos(3*x)

xgrid = np.linspace(0.0, 5.0, 200)
ygrid = f(xgrid)

plt.clf()
plt.plot(xgrid, ygrid)
plt.grid(True)
plt.savefig("IMG_example_5_2.pdf")