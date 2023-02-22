import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

fig, ax = plt.subplots()

xmin, xmax = -5.0, 5.0
ymin, ymax = -5.0, 5.0
Npoints = 200

x = np.linspace(xmin, xmax, Npoints)
y = np.linspace(ymin, ymax, Npoints)
X, Y = np.meshgrid(x, y)

A1 = -1.0
σ1 = 1.0

A2 = 2.0
σ2 = 1.5

arg1 = ((X - 1)**2 + (Y - 3)**2) / σ1**2
arg2 = ((X + 1)**2 + (Y + 2)**2) / σ2**2
Fxy = A1*np.exp(-arg1) + A2*np.exp(-arg2)

CS = ax.contourf(X, Y, Fxy, cmap=cm.coolwarm)
ax.contour(X, Y, Fxy, colors="black")

ax.set_xlabel("Sumbu x")
ax.set_ylabel("Sumbu y")
ax.set_aspect("equal", "box")
ax.set_title("Nama1 (NIM1), Nama2 (NIM2)")

cbar = fig.colorbar(CS)
cbar.ax.set_ylabel("Nilai fungsi")

plt.savefig("IMG_soal_03_kontur.png", dpi=150)
plt.savefig("IMG_soal_03_kontur.pdf")
plt.show()