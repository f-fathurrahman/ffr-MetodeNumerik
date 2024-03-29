import numpy as np
import matplotlib.pyplot as plt

def my_func01(x):
    #return (x - 1) * (x-2)
    #return (x - 1.0) * (x - 2.0) * (x + 3.0)
    #return np.abs(x - 2)
    return np.sin(x)
    #return -np.exp(2*x - x**2)
    #return x**3



def debug_optim_bisection(f, a, b):
    # method parameters
    delta = 0.0001
    SMALL = 0.0002
    NiterMax = 100

    # Initial values for k, n, d
    k = 0 # iteration count
    n = 0 # number of function evaluations

    d = [b - a] # interval length (error)

    # Proposed optimum value
    xOpt = [ (a + b)/2 ] # suboptimal solution

    # Plot
    xgrid = np.linspace(a, b, 200)
    plt.clf()
    plt.plot(xgrid, my_func01(xgrid))
    plt.grid(True)
    plt.title("Optimization with bisection")
    plt.savefig("IMG_func_kiter_{:03d}.png".format(k), dpi=150)

    while (d[-1] >= SMALL) and (k <= NiterMax): # termination criterion

        x1 = (a + b - delta)/2
        x2 = (a + b + delta)/2
    
        # Two function evaluations per iteration
        f1 = my_func01(x1)
        f2 = my_func01(x2)

        #plt.plot([a, b], [my_func01(a), my_func01(b)], marker="o", color="green", linewidth=0)
        #plt.plot([x1, x2], [f1, f2], marker="*", color="red", linewidth=0)

        plt.plot([a, b], [0.0, 0.0], marker="o", color="green", linewidth=0)
        plt.plot([x1, x2], [0.0, 0.0], marker="*", color="red", linewidth=0)

        print()
        print("iteration = ", k)
        print("interval: a = ", a, " b = ", b)
        print("fa = ", my_func01(a), " f(b) = ", my_func01(b))
        print("f1 = ", f1, " f2 = ", f2)

        if f1 <= f2:
            # we have: f(a) > (?)  f(1) <= f(2) < (?) f(b)
            # a is not changed
            b = x2
        else:
            # we have: f(a) > (?)  f(1) > f(2) < (?) f(b)
            # b is not changed
            a = x1
    
        d.append(b - a) # updating interval length
        xOpt.append( (a + b)/2 ) # suboptimal solution

        plt.title("Optimization with bisection")
        plt.savefig("IMG_func_kiter_{:03d}.png".format(k), dpi=150)

        k = k + 1

    if k >= NiterMax:
        print("WARNING: Maximum number of iterations is reached.")



debug_optim_bisection(my_func01, -5.0, 5.0)

