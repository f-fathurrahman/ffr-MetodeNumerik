import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# cm: colormap

x = np.linspace(-5.0, 5.0, 100)
y = np.linspace(-5.0, 5.0, 100)
X, Y = np.meshgrid(x, y)
Z = X + 1j*Y
fZ = np.log(Z)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax and fig will be reused

# First plot
surf1 = ax.plot_surface(X, Y, np.real(fZ), cmap=cm.coolwarm, linewidth=0)
fig.savefig("IMG_meshgrid_01_real.png", dpi=150)

# Second plot, create new fig and plot
ax.cla() # clear axis
surf1 = ax.plot_surface(X, Y, np.imag(fZ), cmap=cm.coolwarm, linewidth=0)
plt.savefig("IMG_meshgrid_01_imag.png", dpi=150)
