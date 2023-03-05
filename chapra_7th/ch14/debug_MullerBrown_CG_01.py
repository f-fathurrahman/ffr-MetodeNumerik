import numpy as np
import matplotlib.pyplot as plt
from muller_brown import *

def my_func(X): # input as vector X
    x, y = X[0], X[1]
    return MullerBrown_energy(x,y)

def my_func_plot(X, Y): # for plotting purpose
    return MullerBrown_energy(X, Y)

def grad_my_func(X): # input as vector X
    x, y = X[0], X[1]
    dfdx, dfdy = MullerBrown_forces(x, y)
    return -np.array([dfdx, dfdy]) # return as numpy array


def linmin_grad(grad_func, x, g, d, αt=3e-5):
    xt = x + αt*d
    gt = grad_func(xt)
    denum = np.dot(g - gt, d)
    if denum != 0.0:
        α = abs( αt * np.dot(g, d)/denum )
    else:
        α = 0.0
    return α


xgrid = np.linspace(-1.5, 1.1, 100)
ygrid = np.linspace(0.0, 2.0, 100)
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)

fig, ax = plt.subplots()
ax.contour(Xgrid, Ygrid, my_func_plot(Xgrid, Ygrid).clip(max=2.0), levels=20)
# clip(max=200) for nicer visualization (isolines)

# Initial point
#x0 = np.array([-1.0, 1.5])
#x0 = np.array([0.0, 1.0])
x0 = np.array([-1.0, 0.75])
NiterMax = 40
func = my_func
grad_func = grad_my_func
x = np.copy(x0)
d_prev = np.zeros(np.size(x0))
g_prev = np.zeros(np.size(x0))

prefix_plt = "IMG_MullerBrown_CG_"
ax.plot(x[0], x[1], marker="o", color="black")
ax.set_aspect("equal")
plt.title("Energy = " + str(func(x)))
plt.savefig(prefix_plt + str(0) + ".png", dpi=150)

for iiter in range(1,NiterMax+1):

    print("\nIteration: ", iiter)
    print("Current point: ", x)

    f = func(x)
    g = grad_func(x)

    if iiter > 1:
        #β = np.dot(g, g)/np.dot(g_prev, g_prev)
        β = np.dot(g-g_prev,g) / np.dot(g_prev,g_prev)
    else:
        β = 0.0

    if β < 0:
        β = 0.0

    #β = 0.0 # set to zero to disable conjugate direction
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

    fnew = func(x)
    iterTryMax = 50
    iterTry = 0
    while fnew > f:
        iterTry += 1
        print("WARNING: Energy does not decrease fnew = ", fnew, " fold = ", f)
        x -= α*d # undo the step
        α *= 0.5 # reduce α
        print("New α = ", α)
        x += α*d
        fnew = func(x)
        if iterTry >= iterTryMax:
            break

    # draw a line from xprev to x
    ax.plot([xprev[0], x[0]], [xprev[1], x[1]], marker="o", color="black")
    plt.title("Energy = " + str(func(x)))
    plt.savefig(prefix_plt + str(iiter) + ".png", dpi=150)

    d_prev = np.copy(d)
    g_prev = np.copy(g)

