# FIXME: Only improper integral is calculated

import sympy
x = sympy.symbols("x")
func_symb = 1/sympy.sqrt(2*sympy.pi)*sympy.exp(-x**2/2)
resExact = sympy.N(sympy.integrate(func_symb, (x, -sympy.oo, -2)))

from integ_routines import *
from math import sqrt, pi, exp

def my_func(x):
    return 1/sqrt(2*pi)*exp(-x**2/2)

a = -100.0 # practical infinity
b =  -2

resApproxInf = sympy.N(sympy.integrate(func_symb, (x, a, b)))

print("\nUsing Boole's rule") # naive
for n in [1, 2, 4, 10, 50, 100]:
    resN = integ_boole_multiple(my_func, a, b, n)
    print("%3d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))


from integ_routines_open import *

def my_func2(t):
    x = 1/t
    return (1/t**2)*my_func(x)

# x -> -oo, t -> 1/x = 0
# x -> -2, t -> 1/x = -1/2
t1 = -0.5
t2 =  0.0

resN = apply_quadrature_multi_interval(
    integ_newtoncotes_open6seg, my_func2, t1, t2, 10
)

print()
print("resN (using Newton-Cotes 6seg, 10 interval)  = ", resN)
print("error (using Newton-Cotes 6seg, 10 interval) = ", abs(resN-resExact))

print()
print("resApproxInf = ", resApproxInf)
print("resExact = ", resExact)
