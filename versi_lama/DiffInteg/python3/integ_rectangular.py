def integ_rectangular(f, x):
    N = len(f)
    assert N == len(x)
    a = x[0]
    b = x[-1]
    s = 0.0
    for i in range(N):
        s = s + f[i]
    return s*(b-a)/N
