import sympy
x = sympy.symbols("x")
func_symb = (x + 1/x)**2

a = 1.0
b = 2.0
resExact = sympy.N(sympy.integrate(func_symb, (x, a, b)))

from integ_romberg import *

def my_func(x):
    return (x + 1/x)**2

resN = integ_romberg(my_func, a, b, MAXIT=8, es=0.5)
print("resN = %18.10f" % resN)
print("res  = %18.10f" % resExact)
print("err  = %18.10e" % abs(resExact-resN))
