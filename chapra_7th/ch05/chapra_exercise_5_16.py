from root_bisection import *
from root_regula_falsi import *
from root_ridder import *

Q = 20.0
g = 9.81

def f(y):
    B = 3 + y
    Ac = 3*y + y**2/2
    return 1 - Q*B/(g*Ac**3)

yl = 0.5
yu = 2.5
yr = root_bisection(f, 0.5, 2.5)
yr = root_regula_falsi(f, 0.5, 2.5)
yr = root_ridder(f, 0.5, 2.5)