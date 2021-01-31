import numpy as np

def calc_quad_root_v1(a, b, c):
    D = b**2 - 4*a*c
    x1 = (-b + np.sqrt(D))/(2*a)
    x2 = (-b - np.sqrt(D))/(2*a)
    return x1, x2

def calc_quad_root_v2(a, b, c):
    D = b**2 - 4*a*c
    x1 = -2*c/(b + np.sqrt(D))
    x2 = -2*c/(b - np.sqrt(D))
    return x1, x2

a = 1.0
b = 3000.001
c = 3.0

x1_true = -0.001
x2_true = -3000.0

x1, x2 = calc_quad_root_v1(a, b, c)
print("Using 1st formula: appprox roots: ", x1, " ", x2)

x1, x2 = calc_quad_root_v2(a, b, c)
print("Using 2nd formula: appprox roots: ", x1, " ", x2)

print("True roots: ", x1_true, " ", x2_true)