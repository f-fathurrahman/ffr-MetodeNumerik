import numpy as np
from integ_trapezoid import *

def integ_romberg( f, a, b, tol=1.0e-6 ):

    def richardson(r,k):
        for j in range(k-1,0,-1):
            c = 4.0**(k-j)
            r[j] = (c*r[j+1] - r[j])/(c-1.0)
        return r

    r = np.zeros(21)