from sympy import *

#u_i_n = Symbol("u_i^{n}")
u_i_n = Symbol("u_i^n")
pprint(u_i_n)

import sympy.printing as printing

delta__y_l = symbols("Delta__y_l")
print(printing.latex(delta__y_l))

u_i_np1 = symbols("u_i__n+1")
print(printing.latex(u_i_np1))

u_ip1_np1 = symbols("u_i+1__n+1")
print(printing.latex(u_ip1_np1))

u_ip1j_np1 = symbols("u_{i+1,j}__n+1")
print(printing.latex(u_ip1j_np1))
pprint(u_ip1j_np1)