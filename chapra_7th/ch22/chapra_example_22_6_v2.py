# FIXME: Only improper integral is calculated

import sympy
SIGMA = 100
x = sympy.symbols("x")
func_symb = 1/(SIGMA*sympy.sqrt(2*sympy.pi))*sympy.exp(-x**2/2/SIGMA**2)
resExact = sympy.N(sympy.integrate(func_symb, (x, -sympy.oo, -2)))

from integ_routines import *
from integ_romberg import *
from math import sqrt, pi, exp
from integ_routines_open import *

def my_func(x):
    return 1/(SIGMA*sqrt(2*pi))*exp(-x**2/2/SIGMA**2)

a = -100000.0 # practical infinity
b = -2

resApproxInf = sympy.N(sympy.integrate(func_symb, (x, a, b)))

print("\nUsing Boole's rule") # naive
for n in [1, 2, 4, 10, 50, 100, 200, 400, 500]:
    resN = apply_quadrature_multi_interval(
        integ_boole, my_func, a, b, n
    )    
    print("%3d %18.10f %18.10e %18.10e" % (n,
        resN, abs(resN-resExact), abs(resN-resApproxInf)))

print("\nUsing integ_romberg")
resN = integ_romberg(my_func, a, b, MAXIT=12)
print("resN = %18.10f %18.10e %18.10e" % (resN,
    abs(resN-resExact), abs(resN-resApproxInf)))

def my_func2(t):
    x = 1/t
    return (1/t**2)*my_func(x)

# x -> -oo, t -> 1/x = 0
# x -> -2, t -> 1/x = -1/2
t1 = 1/b
t2 = 0.0

print("\nUsing Newton-Cotes 6seg")
for n in [1, 2, 4, 10, 50, 100, 200, 500, 1000]:
    resN = apply_quadrature_multi_interval(
        integ_newtoncotes_open6seg, my_func2, t1, t2, n
    )
    print("%5d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))

print()
print("resExact     = %18.12f" % resExact) # exact result from SymPy
print("resApproxInf = %18.12f" % resApproxInf) # result with approximate infinity (a=-100)
print("diff         = %18.10e" % abs(resExact-resApproxInf))