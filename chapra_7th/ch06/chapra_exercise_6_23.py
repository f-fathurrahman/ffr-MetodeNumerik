from math import sqrt

def u(x,y):
    return (x - 4)**2 + (y - 4)**2 - 5

def dudx(x,y):
    return 2*(x - 4)

def dudy(x,y):
    return 2*(y - 4)

def v(x,y):
    return x**2 + y**2 - 16

def dvdx(x,y):
    return 2*x

def dvdy(x,y):
    return 2*y

# Guess
x = 10.0
y = 2.0

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

