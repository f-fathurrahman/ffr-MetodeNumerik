import numpy as np

def f(x):
    return np.exp(-x) - x

def df(x):
    return -np.exp(-x) - 1

x = 0.0

for i in range(1,6):
    xnew = x  - f(x)/df(x) # Newton-Raphson formula
    fxnew = f(xnew)
    print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
    x = xnew
