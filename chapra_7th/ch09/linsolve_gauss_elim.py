import numpy as np

def linsolve_gauss_elim(A_, B_):
    
    N, Nrhs = B_.shape

    A = np.matrix.copy(A_)
    B = np.matrix.copy(B_)

    # Eliminasi
    print("\nForward elimination")
    for k in range(0,N-1):
        for i in range(k+1,N):
            if A[i,k] != 0.0:
                α = A[i,k]/A[k,k]
                print("α = ", α)
                #A[i,k+1:N] = A[i,k+1:N] - alpha*A[k,k+1:N]
                A[i,:] = A[i,:] - α*A[k,:]
                B[i,:] = B[i,:] - α*B[k,:]
            print("Matrix segitiga atas = \n", A)
            print("Vektor b = \n", B)

    print("\nBack substitution ...")
    for k in range(N-1,-1,-1):
        B[k,:] = (B[k,:] - np.dot(A[k,k+1:N],B[k+1:N,:]))/A[k,k]
    print("Done")

    return B # return the matrix, not just the slice

