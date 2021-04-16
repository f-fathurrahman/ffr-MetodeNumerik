import matplotlib.pyplot as plt
import numpy as np

def my_func(x):
    return 1/( (1+x**2)*(1+0.5*x**2) )

x = np.linspace(0.0, 10.0, 1000)

print(my_func(10.0))
print(my_func(100.0))
print(my_func(1000.0))

plt.clf()
plt.plot(x, my_func(x))
plt.title("Chapra Exercise 22.9 (c)")
plt.grid(True)
plt.savefig("IMG_chapra_exercise_22_9_c.pdf")