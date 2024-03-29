import sympy
x = sympy.symbols("x")
func_symb = 1/(x**3 + 1)
a = 1.0
# SymPy cannot evaluate this integral when we used sympy.oo
# We use a large number instead
resExact = sympy.N(sympy.integrate(func_symb, (x, a, 1e10)))

from integ_routines import *
from integ_romberg import *
from math import sqrt, pi, exp
from integ_routines_open import *

def my_func(x):
    return 1/(x**3 + 1)

b = 1000.0 # practical, approximate infinity

resApproxInf = sympy.N(sympy.integrate(func_symb, (x, a, b)))

print()
print("resExact     = %18.12f" % resExact) # exact result from SymPy
print("resApproxInf = %18.12f" % resApproxInf) # result with approximate infinity (a=-100)
print("diff         = %18.10e" % abs(resExact-resApproxInf))

print("\nUsing Boole's rule") # naive
for n in [1, 10, 50, 100, 200, 500, 1000, 2000, 3000, 5000, 10000]:
    resN = apply_quadrature_multi_interval(
        integ_boole, my_func, a, b, n
    )    
    print("%5d %18.10f %18.10e %18.10e" % (n,
        resN, abs(resN-resExact), abs(resN-resApproxInf)))

print("\nUsing integ_romberg")
resN = integ_romberg(my_func, a, b, MAXIT=14, es=1e-12)
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
for n in [1, 2, 4, 10, 50, 100]:
    resN = apply_quadrature_multi_interval(
        integ_newtoncotes_open6seg, my_func2, t1, t2, n
    )
    print("%5d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))

