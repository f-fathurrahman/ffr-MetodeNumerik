from sympy import *

x = symbols("x")
f = exp(-x) - x

print(diff(f, x))