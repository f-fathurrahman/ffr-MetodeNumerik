from math import cos, pi

import sympy
x = sympy.symbols("x")
func_symb = 6 + 3*sympy.cos(x)
resExact = sympy.N(sympy.integrate(func_symb, (x, 0, sympy.pi/2)))

from integ_routines_open import *
from integ_routines import *

def my_func(x):
    return 6 + 3*cos(x)

a = 0.0
b = pi/2

resN = apply_quadrature_multi_interval(
    integ_newtoncotes_open6seg, my_func, a, b, 10
)
print("\nUsing integ_newtoncotes_open6seg 10 intervals")
print("resN = %18.10f" % resN)
print("res  = %18.10f" % resExact)
print("err  = %18.10e" % abs(resExact-resN))

resN = apply_quadrature_multi_interval(
    integ_simpson38, my_func, a, b, 10
)
print("\nUsing integ_simpson38 10 intervals")
print("resN = %18.10f" % resN)
print("res  = %18.10f" % resExact)
print("err  = %18.10e" % abs(resExact-resN))

resN = apply_quadrature_multi_interval(
    integ_boole, my_func, a, b, 10
)
print("\nUsing integ_boole 10 intervals")
print("resN = %18.10f" % resN)
print("res  = %18.10f" % resExact)
print("err  = %18.10e" % abs(resExact-resN))
