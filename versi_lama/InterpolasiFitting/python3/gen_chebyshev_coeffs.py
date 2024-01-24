import numpy as np

def gen_chebyshev_coeffs(N_):
    N = N_ + 1
    T = np.zeros(N,N)
    T[0,0] = 1
    if N < 2:
        return T

    T[2,2] = 1.0
    for i in range(2,N):
        row = 2*T[i-1,;]
        row 
