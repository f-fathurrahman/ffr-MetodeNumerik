import numpy as np
from root_regula_falsi import *

osf1 = 14.621
T1 = 0.0

osf2 = 6.413
T2 = 40.0

def calc_osf(T):
    Ta = T + 273.15
    arg = -139.34411 + 1.575701e5/Ta - 6.642308e7/Ta**2 + 1.2438e10/Ta**3 - 8.621949e11/Ta**4
    return np.exp(arg)

def solve_for_T(osf):
    def f(T):
        Ta = T + 273.15
        return -np.log(osf) - 139.34411 + 1.575701e5/Ta - 6.642308e7/Ta**2 + \
        1.2438e10/Ta**3 - 8.621949e11/Ta**4
    # Initial guess
    T1 = 0.0
    T2 = 40.0
    return root_regula_falsi(f, T1, T2)

Npoints = 10
osf_vals = np.linspace(7.0, 14.0, Npoints)
T = np.zeros(Npoints)
for i,o in enumerate(osf_vals):
    T[i] = solve_for_T(o)

import matplotlib.pyplot as plt
plt.plot(osf_vals, T, marker="o")
plt.grid(True)
plt.xlabel("osf")
plt.ylabel("T")
plt.savefig("IMG_exe_5_18_funcplot.pdf")