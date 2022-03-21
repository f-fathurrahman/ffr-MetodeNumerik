import numpy as np

def newton_interp(x, y, xi):

    assert(len(x) == len(y))
    N = len(x) - 1  # length of array is (N + 1)

    yint = np.zeros(N+1)
    ea = np.zeros(N)

    # finite divided difference table
    fdd = np.zeros((N+1,N+1))

    for i in range(0,N+1):
        fdd[i,0] = y[i]

    for j in range(1,N+1):
        for i in range(0,N-j+1):
            fdd[i,j] = ( fdd[i+1,j-1] - fdd[i,j-1] ) / ( x[i+j] - x[i] )

    xterm = 1.0
    yint[0] = fdd[0,0]
    
    for order in range(1,N+1):
        xterm = xterm * ( xi - x[order-1] )
        yint2 = yint[order-1] + fdd[0,order]*xterm
        ea[order-1] = yint2 - yint[order-1]
        yint[order] = yint2 
    
    return yint, ea
