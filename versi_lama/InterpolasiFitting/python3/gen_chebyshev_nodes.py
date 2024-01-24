import numpy as np

def gen_chebyshev_nodes(N):
    t = np.zeros(N+1)
    for k in range(N+1):
        t[k] = -np.cos(k*np.pi/N)
    return t
