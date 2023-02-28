import numpy as np
from optim_golden_ratio import *

def my_func(x):
    return 2*np.sin(x) - x**2/10

def m_my_func(x):
    return -my_func(x)

xopt, fx = optim_golden_ratio(m_my_func, 0.0, 4.0, TOL=1e-10, verbose=True)

print("\nOptimization result")
print("xopt    = %18.10f" % xopt)
print("f(xopt) = %18.10f" % my_func(xopt))