from math import exp

import sympy
x = sympy.symbols("x")
func_symb = x*sympy.exp(2*x)

a = 0.0
b = 3.0
resExact = sympy.N(sympy.integrate(func_symb, (x, a, b)))

from integ_romberg import *

def my_func(x):
    return x*exp(2*x)

resN = integ_romberg(my_func, a, b, MAXIT=8)
print("resN = %18.10f" % resN)
print("res  = %18.10f" % resExact)
print("err  = %18.10e" % abs(resExact-resN))
