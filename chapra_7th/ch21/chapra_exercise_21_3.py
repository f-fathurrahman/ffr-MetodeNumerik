import sympy
x = sympy.symbols("x")
func_symb = 1 - x - 4*x**3 + 2*x**5
resExact = sympy.N(sympy.integrate(func_symb, (x, -2, 4)))

from integ_routines import *

def my_func(x):
    return 1 - x - 4*x**3 + 2*x**5

a = -2.0
b = 4.0

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
for n in [2, 4]:
    resN = integ_boole_multiple(my_func, a, b, n)
    print("%3d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))
