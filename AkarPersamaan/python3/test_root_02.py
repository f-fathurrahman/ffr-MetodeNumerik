from bisection import *
from regula_falsi import *
from ridder import *
from brent import *

def f(x):
    return x**3 - 10*x**2 + 5.0

x1 = 0.0
x2 = 1.0

x, err = bisection(f, x1, x2, TOL=1e-10, verbose=True)
print("Final root = %18.10f" % x)
print("True err   = %18.10e" % abs(f(x)))

x, err = regula_falsi(f, x1, x2, TOL=1e-10, verbose=True)
print("Final root = %18.10f" % x)
print("True err   = %18.10e" % abs(f(x)))

x, err = ridder(f, x1, x2, TOL=1e-10, verbose=True)
print("Final root = %18.10f" % x)
print("True err   = %18.10e" % abs(f(x)))

x, err = brent(f, x1, x2, TOL=1e-10, verbose=True)
print("Final root = %18.10f" % x)
print("True err   = %18.10e" % abs(f(x)))
