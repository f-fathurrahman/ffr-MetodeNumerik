def my_func(x):
    return 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5

a = 0.0
b = 0.8
f_a = my_func(a)
f_b = my_func(b)

I_exact = 1.640533 # from the book
I = (b - a)*(f_a + f_b)/2
E_t = I_exact - I
ε_t = E_t/I_exact * 100
print("Integral result = %.6f" % I)
print("True error      = %.6f" % E_t)
print("ε_t             = %.1f%%" % ε_t)

import sympy
x = sympy.symbols("x")
f = 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5
d2f = f.diff(x,2)
#sympy.pprint(d2f)
avg_d2f_xi = sympy.integrate( d2f, (x,a,b) )/(b - a)
#print("avg_d2f_xi = ", avg_d2f_xi)
E_a = -1/12*avg_d2f_xi*(b - a)**3
print("Approx error    = %.6f" % E_a)