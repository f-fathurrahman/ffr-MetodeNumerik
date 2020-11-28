import numpy as np
from math import exp

def my_func(t):
    g = 9.8
    m = 68.1
    c = 12.5
    return g*m/c * (1 - exp(-(c/m)*t) )

def calc_exact_sympy():
    # Calculate "Exact" result using SymPy
    import sympy
    g = 9.8
    m = 68.1
    c = 12.5
    t = sympy.symbols("t")
    f = g*m/c * (1 - sympy.exp(-(c/m)*t) )
    d_sympy = sympy.integrate( f, (t,0,10))
    return d_sympy

I_exact = calc_exact_sympy()

a = 0.0
b = 10.0

a_0 = (b + a)/2
a_1 = (b - a)/2
# using a_0 and a_1 as global var
# mapping from x_d -> x
def mapping_func(x_d):
    return a_0 + a_1*x_d
dx_dx_d = (b - a)/2

for NGaussPoints in range(2,7):
    GAUSS_x, GAUSS_c = np.polynomial.legendre.leggauss(NGaussPoints)
    I = 0.0
    for i in range(NGaussPoints):
        x_d = GAUSS_x[i]
        c = GAUSS_c[i]
        I = I + c * my_func( mapping_func(x_d) ) * dx_dx_d
    print("NGaussPoints = %d, I = %.4f" % (NGaussPoints,I))
