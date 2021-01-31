from math import factorial
import numpy as np

def approx_exp(x, N):
    assert(N >= 0)
    if N == 0:
        return np.float32(1.0), np.float32(1.0)
    s = np.float32(0.0)
    for i in range(N+1):
        #s = s + x**np.float32(i)/np.float32(factorial(i))
        s = s + np.float32(x**i)/np.float32(factorial(i))
    #term = x**np.float32(N)/np.float32(factorial(N)) # calculate again
    term = np.float32(x**N)/np.float32(factorial(N)) # calculate again
    return s, term

#x = np.float32(10.0)
x = np.float32(-10.0)
true_val = np.exp(x)
for i in range(50):
    approx_val, term = approx_exp(x, i)
    print("%3d %18.10f %18.10f" % (i+1, approx_val, term))

print(type(approx_val))
print(type(term))

print("true_val = %18.10f" % true_val)
print("error = %10.5e" % abs(true_val - approx_val))