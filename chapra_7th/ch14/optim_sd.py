import numpy as np
import numpy as np

def my_func(X):
    x, y = X[0], X[1]
    return 2*x*y + 2*x - x**2 - 2*y**2

def grad_my_func(X):
    x, y = X[0], X[1]
    dfdx = 2*y + 2 - 2*x
    dfdy = 2*x - 4*y
    return np.array([dfdx, dfdy]) # return as numpy array

def m_my_func(X):
    return -my_func(X)

def grad_m_my_func(X):
    return -grad_my_func(X)



def optim_SD_simple(func, grad_func, x0, α=0.1, NiterMax=1000):
    
    x = np.copy(x0)

    for iiter in range(1,NiterMax+1):

        print("\nIteration: ", iiter)
        print("Current point: ", x)

        f = func(x)
        g = grad_func(x)
        d = -g # step direction, we search for minimum
    
        norm_g = np.sqrt(np.dot(g,g))
        print("f      = %18.10f" % f)
        print("norm g = %18.10e" % norm_g)
        if norm_g < 1e-10:
            print("Converged")
            break

        # Update x
        x = x + α*d 

    return x, f




def linmin_grad(grad_func, x, g, d, αt=1e-5):
    xt = x + αt*d
    gt = grad_func(xt)
    denum = np.dot(g - gt, d)
    if denum != 0.0:
        α = abs( αt * np.dot(g, d)/denum )
    else:
        α = 0.0
    return α


def optim_SD_linmin(func, grad_func, x0, NiterMax=1000):
    x = np.copy(x0)
    for iiter in range(1,NiterMax+1):
        #
        print("\nIteration: ", iiter)
        print("Current point: ", x)
        #
        f = func(x)
        g = grad_func(x)
        d = -g # step direction, we search for minimum
    
        norm_g = np.sqrt(np.dot(g,g))
        print("f      = %18.10f" % f)
        print("norm g = %18.10e" % norm_g)
        if norm_g < 1e-10:
            print("Converged")
            break

        α = linmin_grad(grad_func, x, g, d)
        x = x + α*d

    return x, f




#xopt, fopt = optim_SD_simple(m_my_func, grad_m_my_func, np.array([-1.0, 1.0]))
xopt, fopt = optim_SD_linmin(m_my_func, grad_m_my_func, np.array([-1.0, 1.0]))
print()
print("Optimum point: ", xopt)
print("Optimum value: ", my_func(xopt))



