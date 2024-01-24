import numpy as np

def create_newton_polynom(x, y):
    Ndata = len(x) # jumlah data
    coefs = np.copy(y)
    for k in range(1,Ndata):
        coefs[k:Ndata] = (coefs[k:Ndata] - coefs[k-1])/(x[k:Ndata] - x[k-1])
    return coefs

def eval_newton_polynom(coefs, x, xo):
    N = len(x) - 1 # derajat polinom
    p = coefs[N]
    for k in range(1,N+1):
        p = coefs[N-k] + (xo - x[N-k])*p
    return p

