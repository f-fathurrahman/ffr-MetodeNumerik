import matplotlib.pyplot as plt
import numpy as np

import matplotlib
matplotlib.style.use("ggplot")

t = np.arange(0.0, 2.0, 0.01)
s1 = 1 + np.sin(2*np.pi*t)
s2 = 1 + np.cos(2*np.pi*t)

plt.plot(t, s1, label="s1")
plt.plot(t, s2, label="s2")
plt.xlabel("time (s)")
plt.ylabel("voltage (mV)")
plt.title("About as simple as it gets, folks")
plt.grid(True)
plt.savefig("IMG_02_sin_cos.png", dpi=150)
plt.show()