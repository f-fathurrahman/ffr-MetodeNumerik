import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0.0, 2.0, 0.01)
s1 = 1 + np.sin(2*np.pi*t)
s2 = 1 + np.cos(2*np.pi*t)

fig, ax = plt.subplots()
ax.plot(t, s1, label="s1")
ax.plot(t, s2, label="s2")
ax.set(
    xlabel="time (s)",
    ylabel="voltage (mV)",
    title="About as simple as it gets, folks")
ax.grid(True)
fig.savefig("IMG_01_sin_cos.png", dpi=150)

plt.show()