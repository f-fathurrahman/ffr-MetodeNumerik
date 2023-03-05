import numpy as np
import matplotlib.pyplot as plt

from muller_brown import *

xgrid = np.linspace(-1.5, 0.5, 100)
ygrid = np.linspace(0.0, 2.0, 100)
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)

fig, ax = plt.subplots()
ax.contour(Xgrid, Ygrid, MullerBrown_energy(Xgrid, Ygrid), levels=20)
ax.set_aspect("equal", "box")
plt.savefig("IMG_MullerBrown.png", dpi=150)
