import numpy as np

def calc_quad_root_v1(a, b, c):
    D = np.float32(b**2) - np.float32(4)*a*c
    x1 = (-b + np.sqrt(D))/(np.float32(2)*a)
    x2 = (-b - np.sqrt(D))/(np.float32(2)*a)
    return x1, x2

def calc_quad_root_v2(a, b, c):
    D = np.float32(b**2) - np.float32(4)*a*c
    x1 = -np.float32(2)*c/(b + np.sqrt(D))
    x2 = -np.float32(2)*c/(b - np.sqrt(D))
    return x1, x2

a = np.float32(1.0)
b = np.float32(5000.002)
c = np.float32(10.0)

# [1/500, 5000]
x1_true = np.float32(1/500)
x2_true = np.float32(5000)

x1, x2 = calc_quad_root_v1(a, b, c)
print("Using 1st formula: appprox roots: ", x1, " ", x2)
print(type(x1), type(x2))

x1, x2 = calc_quad_root_v2(a, b, c)
print(type(x1), type(x2))
print("Using 2nd formula: appprox roots: ", x1, " ", x2)

print("True roots: ", x1_true, " ", x2_true)