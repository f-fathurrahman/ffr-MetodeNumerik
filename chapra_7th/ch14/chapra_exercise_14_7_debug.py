import numpy as np
import matplotlib.pyplot as plt

def my_func(X): # input as vector X
    x, y = X[0], X[1]
    return 4*x + 2*y + x**2 - 2*x**4 + 2*x*y - 3*y**2

def my_func_plot(x, y): # for plotting purpose
    return 4*x + 2*y + x**2 - 2*x**4 + 2*x*y - 3*y**2

def grad_my_func(X): # input as vector X
    x, y = X[0], X[1]
    dfdx = 4.0 + 2*x - 8*x**3 + 2*y
    dfdy = 2.0 + 2*x - 6*y
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


xgrid = np.linspace(-3.0, 3.0, 100)
ygrid = np.linspace(-3.0, 3.0, 100)
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)

fig, ax = plt.subplots()
ax.contour(Xgrid, Ygrid, my_func_plot(Xgrid, Ygrid), levels=10)

# Initial point
x0 = np.array([0.0, 0.0])
NiterMax = 40
func = m_my_func
grad_func = grad_m_my_func
x = np.copy(x0)

ax.plot(x[0], x[1], marker="o", color="black")
ax.set_aspect("equal")
plt.savefig("IMG_optim_SD_linmin_exercise_14_7_" + str(0) + ".png", dpi=150)

for iiter in range(1,NiterMax+1):

    print("\nIteration: ", iiter)
    print("Current point: ", x)

    f = func(x)
    g = grad_func(x)
    d = -g # step direction, we search for minimum

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
    plt.savefig("IMG_optim_SD_linmin_exercise_14_7_" + str(iiter) + ".png", dpi=150)

plt.savefig("IMG_chapra_exercise_14_7.pdf")

print("Maximum value = ", my_func(x))