import numpy as np

def gauss_elim(A_, B_):
    
    N, Nrhs = B_.shape

    A = np.matrix.copy(A_)
    B = np.matrix.copy(B_)

    # Eliminasi
    for k in range(0,N-1):
        for i in range(k+1,N):
            if A[i,k] != 0.0:
                alpha = A[i,k]/A[k,k]
                A[i,k+1:N] = A[i,k+1:N] - alpha*A[k,k+1:N]
                B[i,:] = B[i,:] - alpha*B[k,:]
    
    # Substitusi balik
    for k in range(N-1,-1,-1):
        B[k,:] = (B[k,:] - np.dot(A[k,k+1:N],B[k+1:N,:]))/A[k,k]
    
    return B # return the matrix, not just the slice

