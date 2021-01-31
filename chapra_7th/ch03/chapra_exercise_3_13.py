import numpy as np

a = 40.5 # number that we want to approximate the square root

x = 1.0 # guess
x_old = x

NmaxIter = 100
SMALL = np.finfo(np.float64).eps

for i in range(1,NmaxIter+1):
    x = 0.5*( x + a/x )
    err = abs(x-x_old)
    print("%18.10f %10.5e" % (x, err))
    if err <= SMALL:
        print("Convergence is achieved")
        break
    x_old = x

print("x = ", x)
print("true err = ", np.sqrt(a) - x)