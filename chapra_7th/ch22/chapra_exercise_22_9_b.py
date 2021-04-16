import sympy
x = sympy.symbols("x")  # I will use x instead of y
func_symb = sympy.exp(-x)*sympy.sin(x)**2
resExact = sympy.N(sympy.integrate(func_symb, (x, 0, sympy.oo)))

split_point = 6
resExact1 = sympy.N(sympy.integrate(func_symb, (x, split_point, sympy.oo)))
resExact2 = sympy.N(sympy.integrate(func_symb, (x, 0, split_point)))

from integ_routines_open import *
from math import exp, sin

def my_func(x):
    return exp(-x)*sin(x)**2

def my_func2(t):
    x = 1/t
    return (1/t**2)*my_func(x)

# Split the integral into two interval (0,1) and (1,oo)

# The integral (1,oo)

# x -> split_point, t -> 1/x = 1/split_point. split_point != 0
# x -> oo, t -> 1/x = 0
t1 = 0
t2 = 1/split_point

Ninterval = 30
resN1 = apply_quadrature_multi_interval(
    integ_newtoncotes_open6seg, my_func2, t1, t2, Ninterval
)

print()
print("Using Newton-Cotes, Ninterval = ", Ninterval)
print("Result (res,error) = %18.10f %18.10e" % (resN1, abs(resN1-resExact1)))
print("resExact1     = ", resExact1)

from integ_romberg import *

resN2 = integ_romberg(my_func, 0, split_point, MAXIT=8)
print("Using Romberg for another interval")
print("Result (res,error) = %18.10f %18.10e" % (resN2, abs(resN2-resExact2)))

print()
resN = resN1 + resN2
print("Final result: %18.12f %18.12e" % (resN, abs(resN-resExact)))
print("resExact = ", resExact)
