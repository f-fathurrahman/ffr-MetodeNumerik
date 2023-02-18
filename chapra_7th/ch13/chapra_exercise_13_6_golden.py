import numpy as np
from optim_golden_ratio_max import *

def my_func(x):
    return 4*x - 1.8*x**2 + 1.2*x**3 - 0.3*x**4

xopt, fx = optim_golden_ratio(my_func, -2.0, 4.0)
print("xopt    = ", xopt)
print("f(xopt) = ", fx)