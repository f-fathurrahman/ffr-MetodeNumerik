import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.0, 1.0, 200)

Δt = 0.1
A = 1.0
λ = 0.5
f = 2.0
k = -2*np.pi/λ
ω = 2*np.pi*f

t0 = 0.0

fig, ax = plt.subplots()
# ax and fig will be reused

for i in range(20):
    t = t0 + Δt*i
    y = A*np.sin(k*x - ω*t)
    # First plot
    ax.cla()
    ax.plot(x, y)
    ax.grid(True)
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("x")
    fig.savefig("IMG_sin_" + str(i) + ".png", dpi=150)
