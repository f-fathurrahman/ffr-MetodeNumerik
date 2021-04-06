def my_func(x):
    return 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5

# N is number of segments
# Npoints is N + 1
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

a = 0.0
b = 0.8

I_exact = 1.640533 # from the book

print("Using trapezoidal rule")

I_1 = integ_trapz_multiple(my_func, a, b, 1) # satu segmen
ε_t = (I_exact - I_1)/I_exact * 100
print("I_1  = %.6f ε_t = %2.1f%%" % (I_1, ε_t))

I_2 = integ_trapz_multiple(my_func, a, b, 2) # dua segmen
ε_t = (I_exact - I_2)/I_exact * 100
print("I_2  = %.6f ε_t = %2.1f%%" % (I_2, ε_t))

I_3 = integ_trapz_multiple(my_func, a, b, 3) # 3 segmen
ε_t = (I_exact - I_3)/I_exact * 100
print("I_3  = %.6f ε_t = %2.1f%%" % (I_3, ε_t))

I_4 = integ_trapz_multiple(my_func, a, b, 4) # 4 segmen
ε_t = (I_exact - I_4)/I_exact * 100
print("I_4  = %.6f ε_t = %2.1f%%" % (I_4, ε_t))

print()
print("Using Richardson's extrapolation:")

I_12 = 4/3*I_2 - 1/3*I_1
ε_t = (I_exact - I_12)/I_exact * 100
print("I_12 = %.6f ε_t = %2.1f%%" % (I_12, ε_t))


I_24 = 4/3*I_4 - 1/3*I_2
ε_t = (I_exact - I_24)/I_exact * 100
print("I_24 = %.6f ε_t = %2.1f%%" % (I_24, ε_t))

