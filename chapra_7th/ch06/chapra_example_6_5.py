import numpy as np

def f(x):
    return x**10 - 1

def df(x):
    return 10*x**9

x = 0.5

SMALL = np.finfo(np.float64).eps # this is a rather strict convergence criteria

# Try to play with the number of iterations
for i in range(1,100):
    xnew = x  - f(x)/df(x) # Newton-Raphson formula
    fxnew = f(xnew) # should be zero when xnew is a root
    print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
    if np.abs(fxnew) <= SMALL:
        print("Convergence is achieved")
        break
    x = xnew
