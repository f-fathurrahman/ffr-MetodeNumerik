import numpy as np

def LU_decomp(A_):
    
    Nrow, Ncol = A_.shape
    
    assert Nrow == Ncol

    N = Nrow

    A = np.matrix.copy(A_)

    # Eliminasi Gauss (maju)
    for k in range(0,N-1):
        for i in range(k+1,N):
            if A[i,k] != 0.0:
                alpha = A[i,k]/A[k,k]
                A[i,k+1:N] = A[i,k+1:N] - alpha*A[k,k+1:N]
                A[i,k] = alpha
    
    L = np.matrix( np.tril(A,-1) )
    for i in range(N):
        L[i,i] = 1.0 # konstrain Doolittle
    U = np.matrix( np.triu(A) )
    
    return L, U # kembalikan matriks L dan U


A = np.matrix([
[ 2, -2, 0,  0,  0],
[-2, -5, 6,  0,  0],
[ 0, -6, 16, 12, 0],
[ 0,  0, 12, 39, -6],
[ 0,  0, 0,  -6, 14],
],
dtype=np.float64)
print(A)

L, U = LU_decomp(A)
print(L)
print(U)

print(L*U)
print(L*U - A)
