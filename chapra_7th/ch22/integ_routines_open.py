def integ_newtoncotes_open6seg(f, a, b):
    h = (b-a)/6
    x1 = a + h
    x2 = a + 2*h
    x3 = a + 3*h
    x4 = a + 4*h
    x5 = a + 5*h
    #
    s = 11*f(x1) - 14*f(x2) + 26*f(x3) - 14*f(x4) + 11*f(x5)
    return (b-a)*s/20.0


def apply_quadrature_multi_interval(quad_func, f, a, b, Ninterval):
    s = 0.0
    Δ = (b-a)/Ninterval
    for i in range(Ninterval):
        aa = a + i*Δ
        bb = a + (i+1)*Δ
        s = s + quad_func(f, aa, bb)
    return s
