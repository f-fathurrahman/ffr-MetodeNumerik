# ref: Kiusalaas
def integ_trapezeoid_recursive( f, a, b, Iold, k ):
    if k == 1:
        Inew = (f(a) - f(b))*(b - a)/2.0
    else:
        # number of new points
        N = 2**(k-2)
        # 
        h = (b-a)/N
        x = a + h/2.0
        s = 0.0
        for i in range(N):
            s = s + f(x)
            x = x + h
        Inew = (Iold + h*s)/2.0
    return Inew
