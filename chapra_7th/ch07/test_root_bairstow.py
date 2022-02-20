import numpy as np
from numpy.polynomial import Polynomial
from root_bairstow import *

coefs_list = [
    np.array([1.25, -3.875, 2.125, 2.75, -3.5, 1.0]),
    np.array([-2.0, 2.0, -1.0, 1.0]),
    np.array([5.0, -2.0, 6.0, -2.0, 1.0]),
    np.array([8.0, 0.0, 6.0, 0.0, 2.0]),
    np.array([-2, 6.2, -4, 0.7]),
    np.array([9.34, -21.97, 16.3, -3.704]),
    np.array([5.0, -2.0, 6.0, -2.0, 0.0, 10])
]


for coef in coefs_list:
    print("\ncoef = ", coef)
    p = Polynomial(coef)
    print(p)
    print("Using root_bairstow: ", root_bairstow(coef))
    print("Using Numpy:", p.roots())

