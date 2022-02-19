import numpy as np
from numpy.polynomial.polynomial import polydiv

# np.polydiv also can be used, but it is deprecated
a = np.array([-40.0, 2.0, 1.0]) # x^2 + 2*x - 24
d = np.array([-4.0, 1.0]) # x - 4

q, r = polydiv(a, d)
print("a = ", a)
print("d = ", d)
print("q = ", q)
print("r = ", r)
