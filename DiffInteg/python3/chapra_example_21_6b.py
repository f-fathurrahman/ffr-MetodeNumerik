import numpy as np

def my_func(x):
    return 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5

def integ_simpson13( f, a, b ):
    #
    h = (b - a)/2
    x0 = a
    x1 = a + h
    x2 = b
    #
    I = h/3 * ( f(x0) + 4*f(x1) + f(x2) )
    return I

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

import numpy as np
a = 0.0
b = 0.8
Nsegments = 5
x = np.linspace(a, b, Nsegments+1)
I_exact = 1.640533 # from the book

print(x)

I_13 = integ_simpson13(my_func, x[0], x[2])
print("I_13 = %.7f" % I_13 )

I_38 = integ_simpson38(my_func, x[2], x[Nsegments])
print("I_38 = %.7f" % I_38 )

I = I_13 + I_38
E_t = (I_exact - I)
ε_t = E_t/I_exact * 100
print("Integral result = %.7f" % I)
print("True error      = %.7f" % E_t)
print("ε_t             = %.2f%%" % ε_t)
