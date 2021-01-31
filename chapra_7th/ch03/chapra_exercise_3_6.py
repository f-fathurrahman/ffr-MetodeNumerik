from math import factorial
import numpy as np

def calc_exp_minus_v1(x, N):
    assert(x > 0)
    ss = 0.0
    for i in range(N+1):
        ff = (-1)**i
        ss = ss + ff * x**i/factorial(i)
    return ss

def calc_exp_minus_v2(x, N):
    assert(x > 0)
    ss = 0.0
    for i in range(N+1):
        ss = ss + x**i/factorial(i)
    return 1.0/ss

x = 5.0
Nterm = 20

true_val = np.exp(-x)
approx_val_v1 = calc_exp_minus_v1(x, Nterm-1)
approx_val_v2 = calc_exp_minus_v2(x, Nterm-1)

print("err v1 = ", abs(true_val - approx_val_v1))
print("err v2 = ", abs(true_val - approx_val_v2))

# TODO: also calculate relative error