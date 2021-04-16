import matplotlib.pyplot as plt
import numpy as np

def my_func(x):
    return x*np.exp(-x)

x = np.linspace(-2.0, 10.0, 1000)

print(my_func(10.0))
print(my_func(100.0))
print(my_func(1000.0))

plt.clf()
plt.plot(x, my_func(x))
plt.title("Chapra Exercise 22.9 (d)")
plt.grid(True)
plt.savefig("IMG_chapra_exercise_22_9_d.pdf")