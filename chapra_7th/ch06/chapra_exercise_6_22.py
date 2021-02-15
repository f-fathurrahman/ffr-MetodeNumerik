from math import sqrt

def f(x,y):
    u = -x**2 + x + 0.75 - y
    v = y + 5*x*y - x**2
    return u, v

# x <- gx(x,y)
def gx(x,y):
    return y + x**2 - 0.75

# y <- gy(x,y)
def gy(x,y):
    return x**2 - 5*x*y

def gx_v2(x,y):
    return (y + 5*x*y)/x

def gy_v2(x,y):
    return -x**2 + x + 0.75

# Guess root
x = 1.2
y = 1.2

for i in range(1,10):
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