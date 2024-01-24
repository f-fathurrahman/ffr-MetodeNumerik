def integ_trapezoid(f, x):
    N = len(f)
    assert N == len(x)
    a = x[0]
    b = x[-1]
    s = f[0] + f[-1]
    for i in range(1,N-1):
        s = s + 2*f[i]
    return 0.5*s*(b-a)/(N-1)
