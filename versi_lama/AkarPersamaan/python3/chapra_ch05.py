def bisect(f, xl, xu, NiterMax=100, TOL=1e-5):
    
    Niter = 0
    err = 1000.0
    
    xr = 0.5*(xl + xu)
    xr_old = xr

    while True:
                
        Niter = Niter + 1

        tst = f(xl)*f(xr)
        
        if tst < 0.0:
            xu = xr
        elif tst > 0:
            xl = xr
        else:
            err = 0.0

        if (err < TOL) or (Niter > NiterMax):
            print("Breaking the loop: err = ", err)
            break

        xr_old = xr
        xr = 0.5*(xl + xu)

        if xr != 0.0:
            err = abs((xr-xr_old)/xr)

    return xr, err, Niter



