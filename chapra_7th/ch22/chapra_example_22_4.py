from math import sqrt

def my_func(x):
    return 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5

a = 0.0
b = 0.8

a_0 = (b + a)/2
a_1 = (b - a)/2

# using a_0 and a_1 as global var
# mapping from x_d -> x
def mapping_func(x_d):
    return a_0 + a_1*x_d

dx_dx_d = (b - a)/2

# Can be found for example at:
# https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature
NGaussPoints = 3
GAUSS3_c = [5/9, 8/9, 5/9]
GAUSS3_x = [-sqrt(3/5), 0, sqrt(3/5)]

I_exact = 1.640533 # from the book

I = 0.0
for i in range(NGaussPoints):
    x_d = GAUSS3_x[i]
    c = GAUSS3_c[i]
    I = I + c * my_func( mapping_func(x_d) ) * dx_dx_d

ε_t = (I_exact - I)/I_exact * 100
print("I = %.6f, ε_t = %.1f%%" % (I, ε_t))

