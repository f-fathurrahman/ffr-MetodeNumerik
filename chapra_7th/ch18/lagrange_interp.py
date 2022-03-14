def lagrange_interp(x, y, xx):
    assert(len(x) == len(y))
    N = len(x) - 1  # length of array is (N + 1)
    ss = 0.0
    for i in range(0,N+1):
        pp = y[i]
        for j in range(0,N+1):
            if i != j:
                pp = pp*( xx - x[j] ) / ( x[i] - x[j] )
        ss = ss + pp
    return ss
