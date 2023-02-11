import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0.0, 2.2, 0.1)
y = np.arange(0.0, 2.2, 0.1)

X, Y = np.meshgrid(x, y)
u = np.cos(X)*Y
v = np.sin(Y)*Y

fig, ax = plt.subplots()

ax.quiver(X, Y, u, v)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.grid(True)
ax.axis([-0.3, 2.3, -0.3, 2.3])
ax.set_aspect("equal")

plt.show()