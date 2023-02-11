import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-2.0, 2.0, 50)
y = np.linspace(-2.0, 2.0, 50)

X, Y = np.meshgrid(x, y)
Z = X * np.exp(-X**2 - Y**2)

dx, dy = np.gradient(Z) # calculate the gradient numerically

fig, ax = plt.subplots()
ax.contourf(X, Y, Z, levels=10)
ax.quiver(X, Y, dx, dy)

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_aspect("equal", "box")

# show plot
plt.show()
