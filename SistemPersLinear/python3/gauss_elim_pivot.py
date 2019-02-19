import numpy as np

def swap_rows(v,i,j):
    if len(v.shape) == 1:
        # untuk vektor kolom
        v[i],v[j] = v[j],v[i]
    else:
        v[[i,j],:] = v[[j,i],:]

def swap_cols(v,i,j):
    v[:,[i,j]] = v[:,[j,i]]


def gauss_elim_pivot( A_, B_, TOL=1.0e-12 ):
    
    N, Nrhs = B_.shape

    assert Nrhs == 1

    # Bekerja dengan copy agar matriks awal tidak dimodifikasi
    A = np.matrix.copy(A_)
    b = np.matrix.copy(B_[:,0])

    # Faktor skala, nilai absolut terbesar perbaris
    s = np.zeros(N)
    for i in range(N):
        s[i] = np.max(np.abs(A[i,:]))
    
    # ubah s menjadi vektor colom
    s = np.matrix(s).transpose()
    
    for k in range(0,N-1):
        
        # Ganti baris jika diperlukan
        rr = np.abs(A[k:N,k])/s[k:N]
        p = np.argmax(rr) + k
        
        if abs(A[p,k]) < TOL:
            raise RuntimeError("Matriks singular")
        
        if p != k:
            swap_rows(b,k,p)
            swap_rows(s,k,p)
            swap_rows(A,k,p)

        # Eliminasi
        for i in range(k+1,N):
            if A[i,k] != 0.0:
                lam = A[i,k]/A[k,k]
                A[i,k+1:N] = A[i,k+1:N] - lam*A[k,k+1:N]
                b[i] = b[i] - lam*b[k]

    if abs(A[N-1,N-1]) < TOL:
        raise RuntimeError("Matriks singular")
    
    # Substitusi balik
    b[N-1] = b[N-1]/A[N-1,N-1]
    for k in range(N-2,-1,-1):
        b[k] = (b[k] - np.dot(A[k,k+1:N],b[k+1:N]))/A[k,k]
    
    return b

