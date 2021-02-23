import numpy as np

def linsolve_gauss_elim(A_, b_):
    #
    A = np.copy(A_)
    b = np.copy(b_)
    #
    N, _ = A.shape
    s = np.zeros(N)
    for i in range(N):
        s[i] = abs(A[i,0])
        for j in range(1,N):
            if abs(A[i,j]) > s[i]:
                s[i] = abs(A[i,j])
    print("s  = ", s)
    err = _do_elimination(A, b, s)
    if err != 0:
        print("Error when doing elimination")
    print("After elimination:")
    print("A = \n", A)
    #
    x = _do_back_substitution(A,b)
    print("x = \n", x)
    return x

def _do_elimination(A, b, s):
    SMALL = 1e-15
    N, _ = A.shape
    for k in range(N-1):
        #
        _do_pivoting(A, b, s, k)
        #
        if abs(A[k,k]/s[k]) < SMALL:
            return -1
        #
        for i in range(k+1,N):
            #
            factor = A[i,k]/A[k,k]
            #
            for j in range(N): # do for all rows
                A[i,j] = A[i,j] - factor*A[k,j]
            #
            b[i] = b[i] - factor*b[k]
    #
    if abs(A[N-1,N-1]/s[N-1]) < SMALL:
        return -1
    #
    return 0


def _do_pivoting(A, b, s, k):
    N, _ = A.shape
    p = k
    big = abs(A[k,k])
    for ii in range(k+1,N):
        dummy = abs(A[ii,k]/s[ii])
        if dummy > big:
            big = dummy
            p = ii
    if p != k:
        # Swap row
        for jj in range(k,N):
            dummy = A[p,jj]
            A[p,jj] = A[k,jj]
            A[k,jj] = dummy
        # Also swap the RHS
        dummy = b[p]
        b[p] = b[k]
        b[k] = dummy
        #
        dummy = s[p]
        s[p] = s[k]
        s[k] = dummy
    else:
        print("No pivoting is done")
    #
    return


def _do_back_substitution(A, b):
    N, _ = A.shape
    x = np.zeros(N)
    # The last equation
    x[N-1] = b[N-1]/A[N-1,N-1]
    # The remaining equations
    for i in range(N-2,-1,-1): # from N-2, N-1, ..., 0
        ss = 0.0
        for j in range(i+1,N):
            ss = ss + A[i,j]*x[j]
        x[i] = (b[i] - ss)/A[i,i]
    return x
