from bisection import *
from regula_falsi import *

def f(x):
    return 1/((x-0.3)**2 + 0.01) - 1/((x-0.8)**2 + 0.04)

x1 = 0.0
x2 = 1.0

x, err = bisection(f, x1, x2, TOL=1e-10, verbose=True)

print("Final root = %18.10f" % x)
print("err        = %18.10e" % err)

x, err = regula_falsi(f, x1, x2, TOL=1e-10, verbose=True)

print("Final root = %18.10f" % x)
print("err        = %18.10e" % err)
