import numpy as np

def LU_decomp(A_):
    
    N, M = A_.shape
    assert N == M

    A = np.matrix.copy(A_)

    for k in range(N-1):
        for i in range(k+1,N):
            if A[i,k] != 0.0:
                alpha = A[i,k]/A[k,k]
                A[i,k+1:N] = A[i,k+1:N] - alpha*A[k,k+1:N]
                A[i,k] = alpha
    
    return A


def LU_solve(A_, b_):
    
    N, M = A_.shape
    assert N == M

    A = np.matrix.copy(A_)
    b = np.matrix.copy(b_)
    
    for k in range(1,N):
        b[k] = b[k] - np.dot(A[k,0:k],b[0:k])
    b[N-1] = b[N-1]/A[N-1,N-1]
    
    for k in range(N-2,-1,-1):
        b[k] = (b[k] - np.dot(A[k,k+1:N],b[k+1:N]))/A[k,k]
    
    return b
