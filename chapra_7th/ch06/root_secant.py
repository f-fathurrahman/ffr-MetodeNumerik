def root_secant(f, x0, x1, NiterMax=100, TOL=1e-10):
    for i in range(1,NiterMax+1):
        # approximation of derivative of f(x)
        dfx = (f(x0) - f(x1))/(x0 - x1)
        if abs(dfx) < 1e-12:
            print("WARNING: small dfx = ", dfx)
        xnew = x1 - f(x1)/dfx
        #
        fxnew = f(xnew)
        print("secant: %3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(f(xnew)) < TOL:
            x1 = xnew
            break
        x0 = x1
        x1 = xnew
    return x1

def root_secant_mod(f, x1, δ=0.01, NiterMax=100, TOL=1e-10):
    for i in range(1,NiterMax+1):
        # approximation of derivative of f(x)
        dfx = (f(x1+δ) - f(x1))/δ
        if abs(dfx) < 1e-12:
            print("WARNING: small dfx = ", dfx)
        xnew = x1 - f(x1)/dfx
        #
        fxnew = f(xnew)
        print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(f(xnew)) < TOL:
            x1 = xnew
            break
        x1 = xnew
    return x1