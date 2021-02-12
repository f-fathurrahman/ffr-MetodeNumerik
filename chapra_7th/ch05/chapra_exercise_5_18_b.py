import numpy as np
from root_regula_falsi import *

osf1 = 14.621
T1 = 0.0

osf2 = 6.413
T2 = 40.0

def calc_osf(T):
    Ta = T + 273.15
    arg = - 139.34411 + 1.575701e5/Ta - 6.642308e7/Ta**2 + 1.2438e10/Ta**3 - 8.621949e11/Ta**4
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

osf_vals = [8.0, 10.0, 12.0]
for o in osf_vals:
    T = solve_for_T(o)
    osf = calc_osf(T)
    print("abs(osf - o) = ", abs(osf - o))
