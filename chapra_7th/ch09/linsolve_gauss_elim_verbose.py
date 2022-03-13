import numpy as np

def linsolve_gauss_elim(A_, b_, verbose=False):
    
    N = b_.shape[0]

    A = np.copy(A_)
    b = np.copy(b_)

    # Eliminasi
    verbose and print("\nForward elimination")
    for k in range(0,N-1):
        for i in range(k+1,N):
            if A[i,k] != 0.0:
                α = A[i,k]/A[k,k]
                verbose and print("\nα = ", α)
                A[i,:] = A[i,:] - α*A[k,:]
                b[i]   = b[i]   - α*b[k]
            verbose and print("Matrix A = \n", A)
            verbose and print("Vector b = \n", b)

    verbose and print("\nBack substitution ...")
    x = np.zeros((N,1)) # Column vector
    # The last equation
    x[N-1] = b[N-1]/A[N-1,N-1] # second index of x is taken to be 0
    # The remaining equations
    for i in range(N-2,-1,-1): # from N-2, N-1, ..., 0
        ss = 0.0
        for j in range(i+1,N):
            ss = ss + A[i,j]*x[j]
        x[i] = (b[i] - ss)/A[i,i]
    verbose and print("Done")
    verbose and print("Solution x = \n", x)
    return x
