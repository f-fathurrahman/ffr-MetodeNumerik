import numpy as np

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
