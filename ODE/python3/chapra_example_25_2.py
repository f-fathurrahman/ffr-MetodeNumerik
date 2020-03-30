# Local truncation error

# In general it depends of x, y, and h.
# In the present case it only depends on x and h

# In the polynomial case it is easy to calculate the derivatives
# For more complicated cases, SymPy can be used to aid the calculations.

from math import factorial

def trunc_err_second(x,y,h):
    return (-6*x**2 + 24*x - 20)*h**2/factorial(2)

def trunc_err_third(x,y,h):
    return (-12*x + 24)*h**3/factorial(3)

def trunc_err_fourth(x,y,h):
    return -12*h**4/factorial(4)

x = 0.0
y = 1.0 # only needed because it is required as the argument
h = 0.5
E_t2 = trunc_err_second(x,y,h)
E_t3 = trunc_err_third(x,y,h)
E_t4 = trunc_err_fourth(x,y,h)

print("E_t2  = %8.5f" % E_t2)
print("E_t3  = %8.5f" % E_t3)
print("E_t4  = %8.5f" % E_t4)
print("Total = %8.5f" % (E_t2 + E_t3 + E_t4))