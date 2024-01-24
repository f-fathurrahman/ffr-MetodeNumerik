from math import pow, sqrt

def muller(f, a, b, c, TOL=1.0e-9, verbose=False, NiterMax=100 ):

    res = 0.0;
  
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
        a1 = (((d2 * pow(h1, 2)) - 
               (d1 * pow(h2, 2))) / 
              ((h1 * h2) * (h1 - h2))); 
        a2 = (((d1 * h2) - (d2 * h1)) / 
              ((h1 * h2) * (h1 - h2))); 
        x = ((-2 * a0) / (a1 + abs(sqrt(a1 * a1 - 4 * a0 * a2))))
        y = ((-2 * a0) / (a1 - abs(sqrt(a1 * a1 - 4 * a0 * a2))))
  
        # Taking the root which is  
        # closer to x2 
        if (x >= y): 
            res = x + c; 
        else: 
            res = y + c; 
  
        print("%3d %18.10f %18.10e" % (i, res, abs(f(res))))
        #if i > 1:
        #    if abs(res - res_old) <= TOL:
        #        break
        if abs(f(res)) <= TOL:
            break
        
        res_old = res
        a = b; 
        b = c; 
        c = res;

    return res 