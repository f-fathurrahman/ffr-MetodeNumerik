def my_func(x):
    return 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5

# N is number of segments
# Npoints is N + 1
# f is function
def integ_trapz_multiple( f, a, b, N ):
    x0 = a
    xN = b
    h = (b - a)/N
    ss = 0.0
    for i in range(1,N):
        xi = x0 + i*h
        ss = ss + f(xi)
    I = h/2 * ( f(x0) + 2*ss + f(xN) )
    return I

Nsegments = 2
a = 0.0
b = 0.8
I_exact = 1.640533 # from the book

I = integ_trapz_multiple( my_func, a, b, Nsegments )

print("Nsegments = ", Nsegments)
print("Integral result = %.6f" % I)
E_t = I_exact - I
print("True integral   = %.6f" % I_exact)
print("True error      = %.6f" % E_t)
ε_t = E_t/I_exact * 100
print("ε_t             = %.1f%%" % ε_t)

import sympy
x = sympy.symbols("x")
f = 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5
d2f = f.diff(x,2)
#sympy.pprint(d2f)
avg_d2f_xi = sympy.integrate( d2f, (x,a,b) )/(b - a)
#print("avg_d2f_xi = ", avg_d2f_xi)
E_a = -1/12*avg_d2f_xi*(b - a)**3/Nsegments**2
print("Approx error    = %.6f" % E_a)
