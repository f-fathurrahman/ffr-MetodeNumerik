import numpy as np
from optim_golden_ratio_max import *

def my_func(x):
    return 4*x - 1.8*x**2 + 1.2*x**3 - 0.3*x**4

def calc_parabolic_x3(x0, f0, x1, f1, x2, f2):
    # Calculate the optimum here
    num = f0*(x1**2 - x2**2) + f1*(x2**2 - x0**2) + f2*(x0**2 - x1**2)
    denum = 2*( f0*(x1 - x2) + f1*(x2 - x0) + f2*(x0 - x1) )
    return num/denum

# Initial guess
x0 = 1.75; f0 = my_func(x0)
x1 = 2.0; f1 = my_func(x1)
x2 = 2.5; f2 = my_func(x2)

x3 = calc_parabolic_x3(x0, f0, x1, f1, x2, f2)
f3 = my_func(x3)
print("x3 = ", x3, " f3 = ", f3)

SMALL = 1e-10
NiterMax = 100

for iiter in range(1,NiterMax+1):
    xopt_old = x3
    fopt_old = f3

    # Sequentially choose the next points
    x0 = x1; f0 = f1
    x1 = x2; f1 = f2
    x2 = x3; f2 = f3

    x3 = calc_parabolic_x3(x0, f0, x1, f1, x2, f2)
    f3 = my_func(x3)
    print("x3 = ", x3, " f3 = ", f3)

    if abs(fopt_old - f3) < SMALL:
        print("Converged")
        break



