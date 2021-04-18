import sympy

import sympy
x = sympy.symbols("x")
func_symb = 1/(x*(x+2))
a = 2
resExact = sympy.N(sympy.integrate(func_symb, (x, a, sympy.oo)))

from integ_routines import *
from integ_romberg import *
from math import sqrt, pi, exp
from integ_routines_open import *

def my_func(x):
    return 1/(x*(x+2))

b = 1000.0 # practical infinity

resApproxInf = sympy.N(sympy.integrate(func_symb, (x, a, b)))

print("\nUsing Boole's rule") # naive
for n in [1, 2, 4, 10, 50, 100, 200]:
    resN = apply_quadrature_multi_interval(
        integ_boole, my_func, a, b, n
    )    
    print("%5d %18.10f %18.10e %18.10e" % (n,
        resN, abs(resN-resExact), abs(resN-resApproxInf)))

print("\nUsing integ_romberg")
resN = integ_romberg(my_func, a, b, MAXIT=12)
print("resN = %18.10f %18.10e %18.10e" % (resN,
    abs(resN-resExact), abs(resN-resApproxInf)))

def my_func2(t):
    x = 1/t
    return (1/t**2)*my_func(x)

# x -> a, t -> 1/x = 1/a
# x -> oo, t -> 1/x = 0
t1 = 0
t2 = 1/a

print("\nUsing Newton-Cotes 6seg")
for n in [1, 2, 4, 10, 50, 100, 200, 500]:
    resN = apply_quadrature_multi_interval(
        integ_newtoncotes_open6seg, my_func2, t1, t2, n
    )
    print("%5d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))

print()
print("resExact     = %18.12f" % resExact) # exact result from SymPy
print("resApproxInf = %18.12f" % resApproxInf) # result with approximate infinity (a=-100)
print("diff         = %18.10e" % abs(resExact-resApproxInf))