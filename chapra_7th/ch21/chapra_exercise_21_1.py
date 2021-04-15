from math import cos, pi
import numpy as np

import sympy
x = sympy.symbols("x")
func_symb = 6 + 3*sympy.cos(x)
resExact = sympy.N(sympy.integrate(func_symb, (x, 0, sympy.pi/2)))

from integ_routines import *

def my_func(x):
    return 6 + 3*cos(x)

a = 0.0
b = pi/2

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
for n in [2, 4, 10, 20]:
    resN = integ_boole_multiple(my_func, a, b, n)
    print("%3d %18.10f %18.10e" % (n, resN, abs(resN-resExact)))
