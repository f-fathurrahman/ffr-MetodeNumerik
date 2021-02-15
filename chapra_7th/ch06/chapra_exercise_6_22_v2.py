def u(x,y):
    return -x**2 + x + 0.75 - y

def dudx(x,y):
    return -2*x + 1

def dudy(x,y):
    return -1

def v(x,y):
    return y + 5*x*y - x**2

def dvdx(x,y):
    return 5*y - 2*x

def dvdy(x,y):
    return 1 + 5*x

# Guess
x = 1.2
y = 1.2

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
    ff = 0.5*(ui**2 + vi**2)
    print("x, y = %18.10f %18.10f, f = %10.5e" % (x, y, ff))
    # TODO: Check convergence
    if abs(ff) < 1e-10:
        print("converged")
        break
    # Update x
    xnew = x - (ui*J22 - vi*J12)/detJ
    ynew = y - (vi*J11 - ui*J21)/detJ
    x = xnew
    y = ynew

