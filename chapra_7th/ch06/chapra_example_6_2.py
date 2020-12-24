import numpy as np
import matplotlib.pyplot as plt

def f1(x):
    return x

def f2(x):
    return np.exp(-x)

x = np.linspace(0.0, 1.0, 100)
y1 = f1(x)
y2 = f2(x)

plt.clf()
plt.plot(x, y1, label="f1")
plt.plot(x, y2, label="f2")
plt.legend()
plt.savefig("IMG_example_02.pdf")
