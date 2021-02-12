import numpy as np
from root_regula_falsi import *

osf1 = 14.621
T1 = 0.0

osf2 = 6.413
T2 = 40.0

def f1(T):
    osf = osf1
    Ta = T + 273.15
    return -np.log(osf) - 139.34411 + 1.575701e5/Ta - 6.642308e7/Ta**2 + \
    1.2438e10/Ta**3 - 8.621949e11/Ta**4

def f2(T):
    osf = osf2
    Ta = T + 273.15
    return -np.log(osf) - 139.34411 + 1.575701e5/Ta - 6.642308e7/Ta**2 + \
    1.2438e10/Ta**3 - 8.621949e11/Ta**4

root_regula_falsi(f1, -10, 10)

root_regula_falsi(f2, 30, 50)

#print("f = ", f(T1))
