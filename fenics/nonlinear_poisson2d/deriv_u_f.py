from sympy import *

def q(u):
    return 1 + u**2

x, y = symbols("x[0], x[1]")
u = 1 + x + 2*y

f = -diff( q(u)*diff(u,x), x ) - diff( q(u)*diff(u,y), y )
f = simplify(f)

u_code = printing.ccode(u)
f_code = printing.ccode(f)
print(u_code)
print(f_code)

