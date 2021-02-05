import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(1e-10, 10.0, 1000)
y = np.log(x)

plt.plot(x, y)
plt.ylim(-5, 10)
plt.grid(True)
plt.savefig("IMG_ln_x.pdf")