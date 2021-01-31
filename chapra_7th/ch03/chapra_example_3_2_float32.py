from math import factorial
import numpy as np

def approx_exp(x, N):
    assert(N >= 0)
    if N == 0:
        return 1
    s = np.float32(0.0)
    for i in range(N+1):
        s = s + np.float32(x**i)/np.float32(factorial(i))
    return s

x = np.float32(0.5)
true_val = np.exp(x) # from np module

n_digit = 10
# Equation 3.7
ε_s_percent = np.float32(0.5)*np.float32(10**(2-n_digit))

prev_approx = np.float32(0.0)
for N in range(50):
    approx_val = approx_exp(x, N)
    ε_t_percent = abs(approx_val - true_val)/true_val * 100
    if N > 0:
        ε_a_percent = abs(approx_val - prev_approx)/approx_val * 100
    else:
        ε_a_percent = float('nan')
    prev_approx = approx_val
    print("%3d %18.10f %10.5f%% %10.5f%%" % (N+1, approx_val, ε_t_percent, ε_a_percent))
    if ε_a_percent < ε_s_percent:
        print("Converged within %d significant digits" % n_digit)
        break

print("true_val   is %18.10f" % true_val)
print("approx_val is %18.10f" % approx_val)

# Make sure that float32 is used
print()
print("type(true_val)   = ", type(true_val))
print("type(approx_val) = ", type(approx_val))