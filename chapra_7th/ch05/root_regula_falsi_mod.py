def root_regula_falsi_mod(f, xl, xu, NiterMax=100, TOL=1e-10):
    fxl = f(xl)
    fxu = f(xu)
    assert xl < xu
    assert fxl*fxu < 0
    #
    il = 0
    iu = 0
    for i in range(1,NiterMax+1):
        xr = xu - fxu*(xl - xu)/(fxl - fxu)
        fxr = f(xr)
        print("regula_falsi_mod: %4d %18.10f %18.10e" % (i, xr, fxr))
        if abs(fxr) < TOL:
            break # convergence achieved
        if fxl*fxr < 0:
            xu = xr
            fxu = fxu
            iu = 0
            il = il + 1
            if il >= 2:
                fxl = 0.5*fxl
        else:
            xl = xr
            fxl = fxr
            il = 0
            iu = iu + 1
            if iu >= 2:
                fxu = 0.5*fxu
    return xr