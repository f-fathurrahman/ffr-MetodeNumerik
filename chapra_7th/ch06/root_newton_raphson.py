# m is root multiplicity
def root_newton_raphson(f, df, x0, NiterMax=100, TOL=1e-10, m=1):
    x = x0
    for i in range(1,NiterMax+1):
        xnew = x  - m*f(x)/df(x) # Newton-Raphson formula
        fxnew = f(xnew)
        print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(f(xnew)) < TOL:
            x = xnew
            break
        x = xnew
    return x

# Modified Newton-Raphson
# Need 2nd derivative info of f
def root_newton_raphson_mod(f, df, d2f, x0, NiterMax=100, TOL=1e-10):
    x = x0
    for i in range(1,NiterMax+1):
        fx = f(x)
        dfx = df(x)
        d2fx = d2f(x)
        # Using equation 6.16
        xnew = x  - fx*dfx/( dfx**2 - fx*d2fx )
        fxnew = f(xnew)
        print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(f(xnew)) < TOL:
            x = xnew
            break
        x = xnew
    return x