from math import factorial
import numpy as np

def approx_exp(x, N):
    assert(N >= 0)
    if N == 0:
        return 1.0, 1.0
    s = 0.0
    for i in range(N+1):
        s = s + x**i/factorial(i)
    term = x**N/factorial(N) # calculate again
    return s, term

#x = 10.0
x = -10.0
true_val = np.exp(x)
for i in range(50):
    approx_val, term = approx_exp(x, i)
    print("%3d %18.10f %18.10f" % (i+1, approx_val, term))

print("true_val = %18.10f" % true_val)
print("error = %10.5e" % abs(true_val - approx_val))