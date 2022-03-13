import numpy as np

def LU_decomp_pivot(A_):
    
    Nrow, Ncol = A_.shape
    
    assert Nrow == Ncol
    
    N = Nrow

    A = np.matrix.copy(A_)
    
    # Faktor skala
    s = np.zeros((N,1))
    for i in range(N):
        s[i] = np.max(np.abs(A[i,:]))

    iperm = np.arange(N)
    SMALL = np.finfo(np.float64).eps
    # Eliminasi Gauss (maju)
    for k in range(0,N-1):        
        r = np.abs(A[k:N,k])/s[k:N]
        p = np.argmax(r) + k
        if abs(A[p,k]) < SMALL:
            raise RuntimeError("Matriks A singular")
        # Tukar baris jika diperlukan
        if p != k:
            print("INFO: tukar baris %d dengan %d" % (p,k))
            tukar_baris(A, k, p)
            tukar_baris(s, k, p)
            tukar_baris(iperm, k, p)
        
        for i in range(k+1,N):
            if A[i,k] != 0.0:
                α = A[i,k]/A[k,k]
                A[i,k+1:N] = A[i,k+1:N] - α*A[k,k+1:N]
                A[i,k] = α
    
    L = np.matrix( np.tril(A,-1) )
    for i in range(N):
        L[i,i] = 1.0 # konstrain Doolittle
    U = np.matrix( np.triu(A) )
    #
    return L, U, iperm # kembalikan matriks L dan U serta vektor permutasi


def tukar_baris(v, i, j):
    if len(v.shape) == 1: # array satu dimensi atau vektor kolom
        v[i], v[j] = v[j], v[i]
    else:
        v[[i,j],:] = v[[j,i],:]


def LU_solve_pivot(L, U, iperm, b_):
    
    N = L.shape[0]
    
    x = np.zeros((N,1))
    y = np.zeros((N,1))
    
    b = np.copy(b_)
    for i in range(N):
        b[i] = b_[iperm[i]]
    
    # Ly = b
    # Substitusi maju
    y[0] = b[0]/L[0,0]
    for k in range(1,N):
        ss = 0.0
        for j in range(k):
            ss = ss + L[k,j]*y[j]
        y[k] = (b[k] - ss)/L[k,k]
    
    # Ux = y
    # Substitusi balik
    x[N-1] = y[N-1]/U[N-1,N-1]
    for k in range(N-2,-1,-1):
        ss = 0.0
        for j in range(k+1,N):
            ss = ss + U[k,j]*x[j]
        x[k] = (y[k] - ss)/U[k,k]
    
    return x