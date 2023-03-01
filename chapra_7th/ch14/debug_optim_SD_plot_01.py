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

xgrid = np.linspace(-2.0, 4.5, 100)
ygrid = np.linspace(-1.0, 3.0, 100)
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)

fig, ax = plt.subplots()
ax.contour(Xgrid, Ygrid, my_func_plot(Xgrid, Ygrid), levels=10)

# Initial point
x0 = np.array([-1.0, 1.0])
NiterMax = 40
α = 0.1
func = m_my_func
grad_func = grad_m_my_func
x = np.copy(x0)

ax.plot(x[0], x[1], marker="o", color="black")
ax.set_aspect("equal")
plt.savefig("IMG_optim_SD_" + str(0) + ".png", dpi=150)

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
    x = x + α*d
    
    #fprev = f
    #f = func(x)
    #while f > fprev:
    #    x -= α*d # undo step
    #    print("α will be reduced")
    #    α *= 0.9 # reduce α
    #    x += α*d
    #    fprev = f
    #    f = func(x)

    # draw a line from xprev to x
    ax.plot([xprev[0], x[0]], [xprev[1], x[1]], marker="o", color="black")
    plt.savefig("IMG_optim_SD_" + str(iiter) + ".png", dpi=150)

plt.savefig("IMG_debug_optim_SD.pdf")