import numpy as np

def my_func(x):
    return 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5

def integ_simpson38( f, a, b ):
    #
    h = (b - a)/3
    x0 = a
    x1 = a + h
    x2 = a + 2*h
    x3 = b
    #
    I = 3*h/8 * ( f(x0) + 3*f(x1) + 3*f(x2) + f(x3) )
    return I

a = 0.0
b = 0.8
I_exact = 1.640533 # from the book

I = integ_simpson38(my_func, a, b)
E_t = (I_exact - I)
ε_t = E_t/I_exact * 100
print("Integral result = %.7f" % I)
print("True error      = %.7f" % E_t)
print("ε_t             = %.2f%%" % ε_t)

import sympy
x = sympy.symbols("x")
f = 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5
d4f = f.diff(x,4)
avg_d4f_xi = sympy.integrate( d4f, (x,a,b) )/(b - a)
E_a = -avg_d4f_xi*(b - a)**5 / 6480
print("Approx error    = %.7f" % E_a)
