import numpy as np
from math import factorial

Nterms = 4
x = np.pi/3
approx_val = 0.0
for n in range(Nterms+1):
    f = (-1)**(n)
    approx_val = approx_val + f*x**(2*n)/factorial(2*n)

true_val = np.cos(x)

print("true_val = ", true_val)
print("approx_val = ", approx_val)