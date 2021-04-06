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

print("Using Simpson13 rule")

I_2 = integ_simpson13_multiple(my_func, a, b, 2) # dua segmen
ε_t = (I_exact - I_2)/I_exact * 100
print("I_2  = %.6f ε_t = %2.1f%%" % (I_2, ε_t))

I_4 = integ_simpson13_multiple(my_func, a, b, 4) # 4 segmen
ε_t = (I_exact - I_4)/I_exact * 100
print("I_4  = %.6f ε_t = %2.1f%%" % (I_4, ε_t))

print()
print("Using Richardson's extrapolation:")


I_24 = 4/3*I_4 - 1/3*I_2
ε_t = (I_exact - I_24)/I_exact * 100
print("I_24 = %.6f ε_t = %2.1f%%" % (I_24, ε_t))

