from math import sqrt, sin, cos

def u(x,y):
    return x**2 + 1 - y

def dudx(x,y):
    return 2*x

def dudy(x,y):
    return -1

def v(x,y):
    return 2*cos(x) - y

def dvdx(x,y):
    return -2*sin(x)

def dvdy(x,y):
    return -1

# Guess
x = 1.0
y = 1.0

for i in range(1,100):
    # Jacobian matrix elements
    J11 = dudx(x,y)
    J12 = dudy(x,y)
    J21 = dvdx(x,y)
    J22 = dvdy(x,y)
    # 
    detJ = J11*J22 - J12*J21
    #
    ui = u(x,y)
    vi = v(x,y)
    # Update x
    xnew = x - (ui*J22 - vi*J12)/detJ
    ynew = y - (vi*J11 - ui*J21)/detJ

    uval = u(xnew,ynew)
    vval = v(xnew,ynew)
    rmse = sqrt(0.5*(uval**2 + vval**2))

    print("x, y = %18.10f %18.10f %10.5e" % (xnew, ynew, rmse))
    if rmse < 1e-10:
        break
    x = xnew
    y = ynew

