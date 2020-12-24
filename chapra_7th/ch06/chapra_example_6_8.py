import numpy as np

def f(x):
    return np.exp(-x) - x

x0 = 1.0
δ = 0.01

for i in range(1,6):
    # approximation of derivative of f(x)
    dfx = (f(x0+δ) - f(x0))/δ
    #
    xnew = x0 - f(x0)/dfx
    fxnew = f(xnew)
    print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
    x0 = xnew
