import numpy as np

def f(x):
    return np.exp(-x) - x

x0 = 0.0
x1 = 1.0

for i in range(1,6):
    # approximation of derivative of f(x)
    dfx = (f(x0) - f(x1))/(x0 - x1)
    #
    xnew = x1 - f(x1)/dfx
    fxnew = f(xnew)
    print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
    x0 = x1
    x1 = xnew
