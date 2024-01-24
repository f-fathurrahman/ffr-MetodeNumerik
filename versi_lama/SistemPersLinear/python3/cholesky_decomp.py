import numpy as np
import math

def cholesky_decomp(A):
    
    N, M = A.shape
    assert N == M

    L = np.matrix(np.zeros( (3,3) ))
    
    L[0,0] = math.sqrt(A[0,0])

    # Perform the Cholesky decomposition
    for i in range(N):
        for k in range(i+1):
            tmp_sum = 0.0
            for j in range(k):
                tmp_sum = tmp_sum + L[i,j]*L[k,j]
            if (i == k):
                L[i,k] = math.sqrt(A[i,i] - tmp_sum)
            else:
                L[i,k] = (1.0 / L[k,k] * (A[i,k] - tmp_sum))

    return L




def cholesky_solve(L, b_):
    
    N, M = L.shape
    assert N == M

    b = np.matrix.copy(b_)

    # Solusi dari [L]{y} = {b}
    for k in range(N):
        ss = np.dot(L[k,0:k], b[0:k])
        b[k] = (b[k] - ss)/L[k,k]
    
    # Solusi dari [L^T]{x} = {y}
    for k in range(N-1,-1,-1):
        ss = L[k+1:N,k].transpose().dot(b[k+1:N])
        b[k] = (b[k] - ss)/L[k,k]
    
    return b
