import numpy as np

def my_func(x):
    return 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5

# Using equation 21.18
def integ_simpson13_multiple( f, a, b, N ):
    
    assert N >= 2

    x0 = a
    xN = b
    h = (b - a)/N
    x = np.linspace(a, b, N+1)
    
    ss_odd = 0.0
    for i in range(1,N,2):
        ss_odd = ss_odd + f(x[i])
    
    ss_even = 0.0
    for i in range(2,N-1,2):
        ss_even = ss_even + f(x[i])

    I = (b - a)/(3*N) * ( f(x0) + 4*ss_odd + 2*ss_even + f(xN) )
    return I

a = 0.0
b = 0.8
I_exact = 1.640533 # from the book
Nsegments = 2

I = integ_simpson13_multiple(my_func, a, b, Nsegments)
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
E_a = -avg_d4f_xi*(b - a)**5 / (180*Nsegments**4)
print("Approx error    = %.7f" % E_a)
