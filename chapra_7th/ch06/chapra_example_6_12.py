def u(x,y):
    return x**2 + x*y - 10

def dudx(x,y):
    return 2*x + y

def dudy(x,y):
    return x

def v(x,y):
    return y + 3*x*y**2 - 57

def dvdx(x,y):
    return 3*y**2

def dvdy(x,y):
    return 1 + 6*x*y

# Guess
x = 1.5
y = 3.5

for i in range(1,6):
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

    print("x, y = %18.10f %18.10f" % (xnew, ynew))
    # TODO: Check convergence
    x = xnew
    y = ynew

