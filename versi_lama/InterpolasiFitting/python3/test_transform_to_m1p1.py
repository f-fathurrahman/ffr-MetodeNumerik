import numpy as np

from chebyshev_poly import *

def test01():
    x = np.linspace(0.0, 10.0, 10)
    print(x)
    print(transform_to_m1p1(x))

test01()

