import numpy as np
from optim_brent import *

def f(x):
    return 2*np.sin(x) - x**2/10

def mf(x):
    return -f(x)

xopt, fx = optim_brent(mf, 0.0, 4.0, TOL=1e-10)
print("xopt    = ", xopt)
print("f(xopt) = ", fx)