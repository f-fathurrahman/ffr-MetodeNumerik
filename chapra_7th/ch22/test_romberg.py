from math import cos, pi

import sympy
x = sympy.symbols("x")
func_symb = 6 + 3*sympy.cos(x)
resExact = sympy.N(sympy.integrate(func_symb, (x, 0, sympy.pi/2)))

from integ_romberg import *

def my_func(x):
    return 6 + 3*cos(x)

a = 0.0
b = pi/2
resN = integ_romberg(my_func, a, b)
print("resN = %18.12f" % resN)
print("res  = %18.12f" % resExact)
print("err  = %18.12e" % abs(resExact-resN))
