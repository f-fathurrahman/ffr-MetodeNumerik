import numpy as np

def interp_nat_cubic_spline( x, y, d2x, xu ):

    assert len(x) == len(y)
    N = len(x) - 1

    flag = False
    i = 1

    while True:
        #    
        #is_in_interval = (xu >= x[i-1]) and (xu <= x(i))
        is_in_interval = (x[i-1] <= xu <= x[i])
        
        if is_in_interval:
            #
            #write(*,*) 'interval = ', i
            #
            c1 = d2x[i-1]/6.0/( x[i] - x[i-1] )
            c2 = d2x[i]/6.0/( x[i] - x[i-1] )
            c3 = y[i-1]/(x[i] - x[i-1]) - d2x[i-1]*(x[i] - x[i-1])/6.0
            c4 = y[i]/(x[i] - x[i-1]) - d2x[i]*(x[i] - x[i-1])/6.0
            #
            t1 = c1*( x[i] - xu )**3
            t2 = c2*( xu - x[i-1] )**3
            t3 = c3*( x[i] - xu )
            t4 = c4*( xu - x[i-1] )
            yu = t1 + t2 + t3 + t4
            #
            t1 = -3.0*c1*( x[i] - xu )**2
            t2 = 3.0*c2*( xu - x[i-1] )**2
            t3 = -c3
            t4 = c4
            dy = t1 + t2 + t3 + t4
            #
            t1 = 6.0*c1*( x[i] - xu )
            t2 = 6.0*c2*( xu - x[i-1] )
            d2y = t1 + t2
            #
            flag = True
            #
        else:
            #
            i = i + 1
    
        if i == (N + 1) or flag:
            break # break out of the loop

    if not flag:
        raise RuntimeError("xu is outside range of spline")
    
    return yu, dy, d2y


def decomp_trid(e, f, g):
    N = len(e)
    assert N == len(f)
    assert N == len(g)
    #
    for k in range(1,N):
        e[k] = e[k]/f[k-1]
        f[k] = f[k] - e[k]*g[k-1]
    return

# should be called after calling decomp_trid
def subs_trid(e, f, g, r):
    N = len(e)
    assert N == len(f)
    assert N == len(g)

    # Forward subs
    for k in range(1,N):
        r[k] = r[k] - e[k]*r[k-1]
  
    # back subs
    x = np.zeros(N)
    x[N-1] = r[N-1]/f[N-1]
    for k in range(N-2,-1,-1):
        x[k] = ( r[k] - g[k]*x[k+1] ) / f[k]
    return x



def gen_trid_matrix(x, y):

    assert len(x) == len(y)
    N = len(x) - 1

    e = np.zeros(N-1)
    f = np.zeros(N-1)
    g = np.zeros(N-1)
    r = np.zeros(N-1)

    f[0] = 2.0*( x[2] - x[0] )
    g[0] = x[2] - x[1]
    r[0] = 6.0/( x[2] - x[1] ) * ( y[2] - y[1] )
    r[0] = r[0] + 6.0/( x[1] - x[0] ) * ( y[0] - y[1] )
    #
    for i in range(2,N-1):  # in Fortran: 2 until N-2
        e[i-1] = x[i] - x[i-1]
        f[i-1] = 2.0*( x[i+1] - x[i-1] )
        g[i-1] = x[i+1] - x[i]
        r[i-1] = 6.0/( x[i+1] - x[i] ) * ( y[i+1] - y[i] )
        r[i-1] = r[i-1] + 6.0/( x[i] - x[i-1] ) * ( y[i-1] - y[i] )
    
    e[N-2] = x[N-1] - x[N-2]
    f[N-2] = 2.0*( x[N] - x[N-2] )
    r[N-2] = 6.0/( x[N] - x[N-1] ) * ( y[N] - y[N-1] )
    r[N-2] = r[N-2] + 6.0/( x[N-1] - x[N-2] ) * ( y[N-2] - y[N-1] )

    return e, f, g, r
