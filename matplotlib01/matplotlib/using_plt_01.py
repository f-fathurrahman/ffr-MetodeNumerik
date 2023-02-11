import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-5.0, 5.0)
y = np.sin(2*x)

plt.clf()
plt.plot(x, y)
plt.show()

plt.clf()
plt.plot(x, y**2)
plt.show()