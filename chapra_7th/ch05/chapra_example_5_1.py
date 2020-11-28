# Graphical approach

import numpy as np
import matplotlib.pyplot as plt

m = 68.1 # mass, kg
v = 40.0 # velocity, m/s
t = 10.0 # time, s
g = 9.81

def f(c):
    return g*m/c * (1 - np.exp(-c/m*t)) - 40.0

x = np.linspace(4, 20, 200)
y = f(x)
plt.clf()
plt.plot(x, y)

for c in [4,8,12,16,20]:
    print("%4.f %10.3f" % (c, f(c)))
    plt.plot([c], [f(c)], marker="o", color="black")

c = 14.75 # root estimation by graphical inspection
print("At c = %f, f = %f" % (c, f(c)))
plt.plot([c], [f(c)], marker="o", color="black")

plt.gca().set_aspect(0.5)
plt.grid(True)
plt.xlim(0.0, 20.0)
plt.savefig("IMG_example_5_1.pdf")
