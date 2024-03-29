import cmath

# We work with complex arithmetic
def root_muller(f, xr, h=0.5, TOL=1.0e-10, NiterMax=100):

    res = 0.0
    a = complex(xr - h)
    b = complex(xr)
    c = complex(xr + h)

    for i in range(1,NiterMax+1): 
      
        # Calculating various constants  
        # required to calculate x3 
        f1 = f(a)
        f2 = f(b)
        f3 = f(c)
        
        d1 = f1 - f3
        d2 = f2 - f3
        
        h1 = a - c
        h2 = b - c
        
        a0 = f3
        a1 = (((d2 * h1**2) - (d1 * h2**2)) / ((h1 * h2) * (h1 - h2))) 
        a2 = (((d1 * h2) - (d2 * h1)) / ( h1*h2*(h1 - h2) ))
        discr = a1 * a1 - 4 * a0 * a2
        x = -2 * a0 / (a1 + cmath.sqrt(discr))
        y = -2 * a0 / (a1 - cmath.sqrt(discr))
        
        # Taking the root which is closer to x2
        # FIXME: need to check this!!!
        dx = x - c
        dy = y - c
        if abs(dx) >= abs(y-c): 
            res = x + c
        else: 
            res = y + c 

  
        print("muller: %3d %18.10e" % (i, abs(f(res))))
        if abs(f(res)) <= TOL:
            break
        
        res_old = res
        a = b 
        b = c 
        c = res

    return res
