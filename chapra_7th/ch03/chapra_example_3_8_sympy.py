from sympy import *

def calc_quad_root_v1(a, b, c):
    D = b**2 - 4*a*c
    x1 = (-b + sqrt(D))/(2*a)
    x2 = (-b - sqrt(D))/(2*a)
    return x1, x2

def calc_quad_root_v2(a, b, c):
    D = b**2 - 4*a*c
    x1 = -2*c/(b + sqrt(D))
    x2 = -2*c/(b - sqrt(D))
    return x1, x2

a = Rational(1)
b = Rational(3000001, 1000)
c = Rational(3)

x1_true = -Rational(1, 1000)
x2_true = -3000

x1, x2 = calc_quad_root_v1(a, b, c)
print("Using 1st formula: appprox roots: ", x1, " ", x2)
print(type(x1), type(x2))

x1, x2 = calc_quad_root_v2(a, b, c)
print("Using 2nd formula: appprox roots: ", x1, " ", x2)
print(type(x1), type(x2))

print("True roots: ", x1_true, " ", x2_true)