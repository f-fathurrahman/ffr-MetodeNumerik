import numpy as np
from optim_SD_linmin import *

def my_func(X):
    x, y = X[0], X[1]
    return -8*x + x**2 + 12*y + 4*y**2 - 2*x*y

def grad_my_func(X):
    x, y = X[0], X[1]
    dfdx = -8 + 2*x - 2*y
    dfdy = 12 + 8*y - 2*x
    return np.array([dfdx, dfdy]) # return as numpy array


xopt, fopt = optim_SD_linmin(my_func, grad_my_func, np.array([0.0, 0.0]))
print()
print("Optimum point: ", xopt)
print("Optimum value: ", my_func(xopt))
