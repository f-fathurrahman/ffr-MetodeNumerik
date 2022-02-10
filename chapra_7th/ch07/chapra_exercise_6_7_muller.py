import numpy as np
from math import sqrt
import matplotlib.pyplot as plt

def f(x):
    return np.sin(x) + np.cos(1 + x**2) - 1

def root_muller(f, x0, h=0.5, TOL=1.0e-9, verbose=False, NiterMax=100):

    res = 0.0

    a = x0 - h
    b = x0
    c = x0 + h
  
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
        a = b
        b = c
        c = res

    return res

x = np.linspace(0, 2.0, 1000)
fx = f(x)
plt.clf()
plt.plot(x, fx)
plt.grid(True)
plt.savefig("IMG_chapra_exercise_6_7.pdf")

#x1 = 1.0 # complex root?
x1 = 2.0
#x1 = 10.0 # complex
xr = root_muller(f, x1, h=0.5)
print("xr = ", xr)

