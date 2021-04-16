import matplotlib.pyplot as plt
import numpy as np

def my_func(x):
    return 1/(x*(x+2))

x = np.linspace(2.0, 20.0, 1000)
plt.clf()
plt.plot(x, my_func(x))
plt.title("Chapra Exercise 22.9 (a)")
plt.grid(True)
plt.savefig("IMG_chapra_exercise_22_9_a.pdf")