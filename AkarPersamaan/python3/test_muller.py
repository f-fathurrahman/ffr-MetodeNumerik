from muller import *

def f(x): 
    return (1 * pow(x, 3) + 2 * x * x +
                           10 * x - 20); 

# Driver Code 
a = 0.0 
b = 1.0 
c = 2.0 
res = muller(f, a, b, c)

