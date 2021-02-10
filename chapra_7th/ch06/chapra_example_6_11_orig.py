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
def g2(x,y):
    g1 = sqrt(10 - x*y)
    g2 = sqrt( (57 - y)/(3*x) )
    return g1, g2

# Guess root
x = 1.5
y = 3.5

for i in range(1,10):
    xnew, ynew = g2(x,y)
    Δx = abs(xnew - x)
    Δy = abs(ynew - y)
    print("First iter:")
    print("xnew = ", xnew)
    print("ynew = ", ynew)
    print(Δx, " ", Δy)
    x = xnew
    y = ynew

