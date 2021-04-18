import matplotlib.pyplot as plt
import numpy as np

SIGMA = 1000.0
def my_func(x):
    return 1/(SIGMA*np.sqrt(2*np.pi))*np.exp(-x**2/2/SIGMA**2)

x = np.linspace(-5.0, 5.0, 1000)

plt.clf()
plt.plot(x, my_func(x))
plt.title("SIGMA = %f" % SIGMA)
plt.grid(True)
plt.savefig("IMG_gaussian.pdf")