def root_regula_falsi(f, xl, xu, NiterMax=100, TOL=1e-10):
    fxl = f(xl)
    fxu = f(xu)
    assert xl < xu
    assert fxl*fxu < 0
    #
    for i in range(1,NiterMax+1):
        xr = xu - fxu*(xl - xu)/(fxl - fxu)
        fxr = f(xr)
        print("%4d %18.10f %18.10e" % (i, xr, fxr))
        if abs(fxr) < TOL:
            break # convergence achieved
        if fxl*fxr < 0:
            xu = xr
            fxu = fxr
        else:
            xl = xr
            fxl = fxr
    return xr