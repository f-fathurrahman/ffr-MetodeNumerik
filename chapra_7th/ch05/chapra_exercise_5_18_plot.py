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

import matplotlib.pyplot as plt

T = np.linspace(0.0, 40.0, 50)
osf = calc_osf(T)
plt.clf()
plt.plot(T, osf)
plt.xlabel("T")
plt.ylabel("osf")
plt.grid(True)
plt.savefig("IMG_exe_5_18_plot.pdf")