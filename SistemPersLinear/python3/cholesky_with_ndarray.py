import numpy as np
import math

def cholesky_decomp(A_):

    N, M = A_.shape
    
    assert N == M

    A = np.copy(A_)

    for k in range(N):
        try:
            A[k,k] = math.sqrt(A[k,k] - np.dot(A[k,0:k],A[k,0:k]))
        except ValueError:
            raise RuntimeException("Matriks tidak definit positif")
    
        for i in range(k+1,N):
            A[i,k] = (A[i,k] - np.dot(A[i,0:k],A[k,0:k]))/A[k,k]

    return np.tril(A)


def cholesky_solve(L,b_):

    N, _ = L.shape

    b = np.copy(b_)
    
    # Solution of [L]{y} = {b}
    for k in range(N):
        b[k] = (b[k] - np.dot(L[k,0:k],b[0:k]))/L[k,k]
    
    # Solution of [L_transpose]{x} = {y}
    for k in range(N-1,-1,-1):
        b[k] = (b[k] - np.dot(L[k+1:N,k],b[k+1:N]))/L[k,k]
    
    return b


import numpy as np

A = np.array([
    [4, -2, 2],
    [-2, 2, -4],
    [2, -4, 11]], dtype=np.float64)
print(A)

b = np.array([1, 3/2, 3], dtype=np.float64)
print(b)

L = cholesky_decomp(A)
print("L after cholesky_decomp:")
print(L)
print("Check L (should be zeros)")
print(np.matmul(L,L.transpose()) - A)

x = cholesky_solve(L, b)

print("x = ")
print(x)

print("Check result (should be zeros)")
print(np.matmul(A,x)-b)
