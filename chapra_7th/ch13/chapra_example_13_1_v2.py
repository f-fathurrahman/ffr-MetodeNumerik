import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 2*np.sin(x) - x**2/10

def plot_my_func(myf, xl, xu, x1, x2, filesave="IMG_f.png", title=None):
    if abs(xu - xl) > 1e-3:
        xgrid = np.linspace(xl, xu, 200)
    else:
        xgrid = np.linspace(xl - 1e-3, xu + 1e-3, 200)
    ygrid = myf(xgrid)
    plt.clf()
    plt.plot(xgrid, ygrid, color="black")
    plt.plot([xl], [myf(xl)], marker="o", color="red", label="xl")
    plt.plot([xu], [myf(xu)], marker="o", color="red", label="xu")
    plt.plot([x1], [myf(x1)], marker="*", color="blue", label="x1")
    plt.plot([x2], [myf(x2)], marker="*", color="magenta", label="x2")
    if title is not None:
        plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filesave, dpi=150)

xl = 0.0
xu = 4.0
xgrid = np.linspace(xl, xu, 200)
ygrid = f(xgrid)

γ = (np.sqrt(5) - 1)/2

d = γ*(xu - xl)
x1 = xl + d
x2 = xu - d
f1 = f(x1)
f2 = f(x2)

plot_my_func(f, xl, xu, x1, x2, filesave="IMG_iter_0.png")

iiter = 0

SMALL = np.finfo(np.float64).resolution
NiterMax = 100
for iiter in range(1,NiterMax+1):
    xint = xu - xl
    if f2 > f1:
        xu = x1 # x1 becomes upper bound
        x1 = x2
        f1 = f2
        d = γ*(xu - xl)
        x2 = xu - d
        f2 = f(x2)
    else:
        xl = x2
        x2 = x1
        f2 = f1
        d = γ*(xu - xl)
        x1 = xl + d
        f1 = f(x1)

    if f1 > f2:
        xopt = x1
        fx = f1
        title = "x1 is xopt"
    else:
        xopt = x2
        fx = f2
        title = "x2 is xopt"

    plot_my_func(f, xl, xu, x1, x2, filesave="IMG_iter_" + str(iiter) + ".png", title=title)
    print("%3d %18.10f %18.10f %18.10e" % (iiter, xopt, fx, d))

    if abs(xopt) > SMALL:
        ea = (1 - γ)*abs(xint/xopt)

    if ea <= 1e-10:
        print("Converged")
        break
