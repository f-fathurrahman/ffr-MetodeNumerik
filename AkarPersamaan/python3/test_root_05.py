from bisection import *
from regula_falsi import *

def f(x):
    return x**3 - 75

x1 = 3.0
x2 = 5.0

x, err = bisection(f, x1, x2, TOL=1e-10, verbose=True)

print("Final root = %18.10f" % x)
print("err        = %18.10e" % err)

x, err = regula_falsi(f, x1, x2, TOL=1e-10, verbose=True)

print("Final root = %18.10f" % x)
print("err        = %18.10e" % err)
