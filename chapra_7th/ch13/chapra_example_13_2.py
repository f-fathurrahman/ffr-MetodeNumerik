import numpy as np

def my_func(x):
    return 2*np.sin(x) - x**2/10

# Given three points and function f, this will fit a parabola
# to f and calculate its optimum value
def calc_parabolic_x3_with_f(f, x0, x1, x2):
    # We evaluate f0, f1, and f2 here
    # In case function evaluation is expensive this can be made into
    # arguments instead
    f0 = f(x0)
    f1 = f(x1)
    f2 = f(x2)
    # Calculate the optimum here
    num = f0*(x1**2 - x2**2) + f1*(x2**2 - x0**2) + f2*(x0**2 - x1**2)
    denum = 2*( f0*(x1 - x2) + f1*(x2 - x0) + f2*(x0 - x1) )
    return num/denum


def calc_parabolic_x3(x0, f0, x1, f1, x2, f2):
    # Calculate the optimum here
    num = f0*(x1**2 - x2**2) + f1*(x2**2 - x0**2) + f2*(x0**2 - x1**2)
    denum = 2*( f0*(x1 - x2) + f1*(x2 - x0) + f2*(x0 - x1) )
    return num/denum

# Initial guess
x0 = 0.0; f0 = my_func(x0)
x1 = 1.0; f1 = my_func(x1)
x2 = 4.0; f2 = my_func(x2)

x3 = calc_parabolic_x3(x0, f0, x1, f1, x2, f2)
f3 = my_func(x3)
print("x3 = %18.10f f3 = %18.10f" % (x3, f3))

TOL = 1e-10
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
    print("x3 = %18.10f f3 = %18.10f" % (x3, f3))

    if abs(fopt_old - f3) < TOL:
        print("Converged")
        break



