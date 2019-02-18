import numpy as np
from gauss_elim import *

def vandermode(v):
    n = len(v)
    a = np.zeros((n,n))
    for j in range(n):
        a[:,j] = v**(n-j-1)
    return a

v = np.array([1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
b = np.matrix([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]).transpose()

A = np.matrix(vandermode(v))
print(A)

x = gauss_elim(A,b)

det = np.prod(np.diagonal(A))

print("x =\n",x)
print("det = \n",det)
print("Check result: [a]{x} - b =\n", np.dot(A,x) - b)
