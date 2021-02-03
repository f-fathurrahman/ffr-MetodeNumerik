import numpy as np

def calc_quad_root_v1(a, b, c):
    D = np.float128(b**2) - np.float128(4)*a*c
    x1 = (-b + np.sqrt(D))/(np.float128(2)*a)
    x2 = (-b - np.sqrt(D))/(np.float128(2)*a)
    return x1, x2

def calc_quad_root_v2(a, b, c):
    D = np.float128(b**2) - np.float128(4)*a*c
    x1 = -np.float128(2)*c/(b + np.sqrt(D))
    x2 = -np.float128(2)*c/(b - np.sqrt(D))
    return x1, x2

a = np.float128(1.0)
b = np.float128(3000.001)
c = np.float128(3.0)

x1_true = np.float128(-0.001)
x2_true = np.float128(-3000.0)

x1, x2 = calc_quad_root_v1(a, b, c)
print("Using 1st formula: approx roots: ", x1, " ", x2)
print(type(x1), type(x2))

x1, x2 = calc_quad_root_v2(a, b, c)
print(type(x1), type(x2))
print("Using 2nd formula: approx roots: ", x1, " ", x2)

print("True roots: ", x1_true, " ", x2_true)