from my_root_finders import root_bisection

def my_func(T):
    cp = 1.2
    a0 = 0.99403
    a1 = 1.671e-4
    a2 = 9.7215e-8
    a3 = -9.5838e-11
    a4 = 1.952e-14
    return a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4 - cp

import matplotlib.pyplot as plt
import numpy as np

Tgrid = np.linspace(900, 1500)
f = my_func(Tgrid)
plt.clf()
plt.plot(Tgrid, f)
plt.grid(True)
plt.savefig("IMG_chapra_exercise_8_36.pdf")

xroot = root_bisection(my_func, 1100.0, 1200.0, verbose=True)