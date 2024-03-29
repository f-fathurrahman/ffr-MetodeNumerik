import numpy as np
import matplotlib.pyplot as plt

def my_func01(x):
    return (x - 1) * (x-2)
    #return (x - 1.0) * (x - 2.0) * (x + 3.0)
    #return np.abs(x - 2)
    #return np.sin(x)
    #return -np.exp(2*x - x**2)
    #return x**3

γ = (-1.0 + np.sqrt(5))/2

# method parameters
delta = 0.0001
SMALL = 0.0002


# Initial values for k, n, d
k = 0 # iteration count
n = 0 # number of function evaluations

a = -5.0
b = 5.0
d = [b - a] # interval length (error)

# optimization loop
xOpt = [ (a + b)/2 ] # suboptimal solution

# Plot
xgrid = np.linspace(a, b, 200)
plt.clf()
plt.plot(xgrid, my_func01(xgrid))
plt.grid(True)
plt.savefig("IMG_func_00.png", dpi=150)

f1 = np.nan
f2 = np.nan

NiterMax = 100
while (d[-1] >= SMALL) and (k <= NiterMax): # termination criterion

    x1 = γ*a + (1 - γ)*b
    x2 = (1 - γ)*a + γ*b
    
    if k == 1:
        f1 = my_func01(x1)
        f2 = my_func01(x2) # two function evaluations
    elif abs(x2 - xOpt[-1]) < abs(x1 - xOpt[-1]):
        f2 = f1
        f1 = my_func01(x1) # one function evaluation
    else:
        f1 = f2
        f2 = my_func01(x2) # one function evaluation

    #plt.plot([a, b], [my_func01(a), my_func01(b)], marker="o", color="green", linewidth=0)
    #plt.plot([x1, x2], [f1, f2], marker="*", color="red", linewidth=0)

    # Only plot the abscissa, the ordinate is set to zero
    plt.plot([a, b], [0.0, 0.0], marker="o", color="green", linewidth=0)
    plt.plot([x1, x2], [0.0, 0.0], marker="*", color="red", linewidth=0)

    plt.savefig("IMG_func_k_" + str(k) + ".png", dpi=150)

    print()
    print("iteration = ", k)
    print("interval: a = ", a, " b = ", b)
    print("fa = ", my_func01(a), " f(b) = ", my_func01(b))
    print("f1 = ", f1, " f2 = ", f2)

    if f1 <= f2:
        b = x2
        xOpt.append(x1)
    else:
        a = x1
        xOpt.append(x2)

    d.append(b-a)

    k = k + 1


if k >= NiterMax:
    print("WARNING: Maximum number of iterations is reached.")
