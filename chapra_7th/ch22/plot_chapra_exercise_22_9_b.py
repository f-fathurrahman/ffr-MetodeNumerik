import matplotlib.pyplot as plt
import numpy as np

def my_func(x):
    return np.exp(-x)*np.sin(x)**2

x = np.linspace(0.0, 10.0, 1000)
plt.clf()
plt.plot(x, my_func(x))
plt.title("Chapra Exercise 22.9 (b)")
plt.grid(True)
plt.savefig("IMG_chapra_exercise_22_9_b.pdf")