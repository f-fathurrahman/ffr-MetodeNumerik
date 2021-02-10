from math import sqrt

def f(x,y):
    u = x**2 + x*y - 10
    v = y + 3*x*y**2 - 57
    return u, v

def g(x,y):
    g1 = (10 - x**2)/y
    g2 = 57 - 3*x*y**2
    return g1, g2

# 2nd formulation
def gx_v2(x,y):
    return sqrt(10 - x*y)

def gy_v2(x,y):
    return sqrt( (57 - y)/(3*x) )

# Guess root
x = 1.5
y = 3.5

for i in range(1,20):
    xnew = gx_v2(x,y)
    Δx = abs(x - xnew)
    x = xnew
    ynew = gy_v2(x,y)
    Δy = abs(y - ynew)
    y = ynew
    print("Iter: ", i)
    print("x, y = %18.10f %18.10f" % (xnew, ynew))
    print("Δx, Δy = %10.5e %10.5e" % (Δx, Δy))

print("At root: ", f(x,y), " (should be close to zero)")