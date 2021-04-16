import sympy
x = sympy.symbols("x")
func_symb = 1/(x*(x+2))
resExact = sympy.N(sympy.integrate(func_symb, (x, 2, sympy.oo)))

from integ_routines_open import *

def my_func(x):
    return 1/(x*(x+2))

def my_func2(t):
    x = 1/t
    return (1/t**2)*my_func(x)

# x -> 2, t -> 1/x = 1/2
# x -> oo, t -> 1/x = 0
t1 = 0
t2 = 1.0/2.0

Ninterval = 30
resN = apply_quadrature_multi_interval(
    integ_newtoncotes_open6seg, my_func2, t1, t2, Ninterval
)

print()
print("Using Newton-Cotes, Ninterval = ", Ninterval)
print("Result (res,error) = %18.10f %18.10e" % (resN, abs(resN-resExact)))
print("resExact     = ", resExact)
