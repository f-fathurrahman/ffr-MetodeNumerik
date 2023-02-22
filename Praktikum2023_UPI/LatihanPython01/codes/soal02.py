import numpy as np
import matplotlib.pyplot as plt

xmin = -5.0
xmax = 5.0
Npoints = 500

x = np.linspace(xmin, xmax, Npoints)
y1 = np.sin(x)
y2 = np.cos(x)

plt.clf()
plt.plot(x, y1, label="sin(x)")
plt.plot(x, y2, label="cos(x)")
plt.legend()
plt.xlabel("Sumbu x")
plt.ylabel("Sumbu y")
plt.title("Nama1 (NIM1), Nama2 (NIM2), Nama3 (NIM3)")
plt.grid(True)
plt.savefig("IMG_sin_cos_01.png", dpi=150)
plt.savefig("IMG_sin_cos_01.jpg", dpi=150)
plt.savefig("IMG_sin_cos_01.pdf")
plt.show()
