import numpy as np
import matplotlib.pyplot as plt

def my_func(X): # input as vector X
    x, y = X[0], X[1]
    return 2*x*y + 2*x - x**2 - 2*y**2

def my_func_plot(X, Y): # for plotting purpose
    return 2*X*Y + 2*X - X**2 - 2*Y**2    

def grad_my_func(X): # input as vector X
    x, y = X[0], X[1]
    dfdx = 2*y + 2 - 2*x
    dfdy = 2*x - 4*y
    return np.array([dfdx, dfdy]) # return as numpy array

def m_my_func(X):
    return -my_func(X)

def grad_m_my_func(X):
    return -grad_my_func(X)

def linmin_grad(grad_func, x, g, d, αt=1e-5):
    xt = x + αt*d
    gt = grad_func(xt)
    denum = np.dot(g - gt, d)
    if denum != 0.0:
        α = abs( αt * np.dot(g, d)/denum )
    else:
        α = 0.0
    return α


xgrid = np.linspace(-2.0, 4.5, 100)
ygrid = np.linspace(-1.0, 3.0, 100)
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)

fig, ax = plt.subplots()
ax.contour(Xgrid, Ygrid, my_func_plot(Xgrid, Ygrid), levels=10)

# Initial point
x0 = np.array([-1.0, 1.0])
NiterMax = 40
func = m_my_func
grad_func = grad_m_my_func
x = np.copy(x0)
d_prev = np.zeros(np.size(x0))
g_prev = np.zeros(np.size(x0))

ax.plot(x[0], x[1], marker="o", color="black")
ax.set_aspect("equal")
plt.savefig("IMG_optim_CG_linmin_" + str(0) + ".png", dpi=150)

for iiter in range(1,NiterMax+1):

    print("\nIteration: ", iiter)
    print("Current point: ", x)

    f = func(x)
    g = grad_func(x)

    if iiter > 1:
        β = np.dot(g, g)/np.dot(g_prev, g_prev)
        #β = np.dot(g-g_prev,g) / np.dot(g_prev,g_prev)
    else:
        β = 0.0

    if β < 0:
        β = 0.0

    print("β = ", β)
    d = -g + β*d_prev # step direction, we search for minimum

    ax.quiver(x[0], x[1], d[0], d[1], color="blue") # also plot the direction

    norm_g = np.sqrt(np.dot(g,g))
    print("f      = %18.10f" % f)
    print("norm g = %18.10e" % norm_g)
    if norm_g < 1e-10:
        print("Converged")
        break

    # Update x
    xprev = np.copy(x)
    α = linmin_grad(grad_func, x, g, d)
    x = x + α*d 
    # draw a line from xprev to x
    ax.plot([xprev[0], x[0]], [xprev[1], x[1]], marker="o", color="black")
    plt.savefig("IMG_optim_CG_linmin_" + str(iiter) + ".png", dpi=150)

    d_prev = np.copy(d)
    g_prev = np.copy(g)

