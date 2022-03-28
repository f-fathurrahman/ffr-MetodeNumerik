import numpy as np
from optim_golden_ratio import *

def f(x):
    return 2*np.sin(x) - x**2/10

xopt, fx = optim_golden_ratio(f, 0.0, 4.0)
print("xopt    = ", xopt)
print("f(xopt) = ", fx)