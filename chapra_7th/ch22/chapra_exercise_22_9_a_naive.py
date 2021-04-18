import sympy
x = sympy.symbols("x")
func_symb = 1/(x*(x+2))
resExact = sympy.N(sympy.integrate(func_symb, (x, 2, sympy.oo)))

from integ_routines import *
from integ_romberg import integ_romberg

def my_func(x):
    return 1/(x*(x+2))

a = 2.0
b = 1000.0 # practical infinity

resApproxInf = sympy.N(sympy.integrate(func_symb, (x, 2, b)))

print("\nUsing trapz")
for n in [1, 2, 4, 10, 50, 100]:
    resN = integ_trapz_multiple(my_func, a, b, n)
    print("%3d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))

print("\nUsing simpson13")
for n in [2, 4, 10, 50, 100]:
    resN = integ_simpson13_multiple(my_func, a, b, n)
    print("%3d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))

print("\nUsing simpson13_v2")
for n in [1,2,10,100]:
    resN = integ_simpson13_multiple_v2(my_func, a, b, n)
    print("%3d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))

print("\nUsing simpson38")
for n in [2, 4, 10, 50, 100]:
    resN = integ_simpson38_multiple(my_func, a, b, n)
    print("%3d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))

print("\nUsing Boole")
for n in [2, 4, 10, 20, 50, 100, 200, 500]:
    resN = integ_boole_multiple(my_func, a, b, n)
    print("%3d %18.10f %18.10e %18.10e" % (n, resN, abs(resN-resExact), abs(resN-resApproxInf)))

print("\nUsing integ_romberg")
resN = integ_romberg(my_func, a, b, MAXIT=15)
print("resN = %18.10f %18.10e" % (resN, abs(resN-resExact)))
print("resN w.r.t approxInf = %18.10f %18.10e" % (resN, abs(resN-resApproxInf)))

print()
print("resExact     = ", resExact)
print("resApproxInf = ", resApproxInf)
print("diff         = ", abs(resExact-resApproxInf))