import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 2*np.sin(x) - x**2/10

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

iiter = 0
plt.clf()
plt.plot(xgrid, ygrid)
plt.plot([xl], [f(xl)], marker="o", label="l")
plt.plot([xu], [f(xu)], marker="o", label="u")
plt.plot([x1], [f1], marker="*", label="1")
plt.plot([x2], [f2], marker="^", label="2")
plt.legend()
plt.tight_layout()
plt.savefig("IMG_chapra_example_13_1_" + str(iiter) + ".png", dpi=150)

SMALL = np.finfo(np.float64).resolution
NiterMax = 100
for iiter in range(NiterMax):
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
    else:
        xopt = x2
        fx = f2
    print("%3d %18.10f %18.10f %18.10e" % (iiter, xopt, fx, d))

    if abs(xopt) > SMALL:
        ea = (1 - γ)*abs(xint/xopt)

    if ea <= 1e-10:
        print("Converged")
        break
