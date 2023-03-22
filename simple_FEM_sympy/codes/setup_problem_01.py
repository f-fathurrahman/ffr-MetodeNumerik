from sympy import *

init_printing(use_unicode=True)

# Spatial variable
x = Symbol("x", real=True)

# Domain length
L = Symbol("L", real=True, positive=True)


# Setup linear element
Nelements = 2 
Nnodes = Nelements + 1

# Setup the grid
xnodes = [Integer(0), L/2, L] # I use Integer(0) to make sure that all xnodes are SymPy objects
h = [xnodes[i+1]-xnodes[i] for i in range(Nnodes-1)]

assert(len(xnodes) == Nnodes)
assert(len(h) == Nelements)


# Setup basis functions
Nfuncs = []

# First
cond1 = (x >= xnodes[0]) & (x <= xnodes[1])
f1 = (xnodes[1] - x)/h[0]
Nfuncs.append( Piecewise( (f1, cond1), (0, True) ) )

# Second
cond1 = (x >= xnodes[0]) & (x <= xnodes[1])
f1 = ( x - xnodes[0])/h[0]
#
cond2 = (x >= xnodes[1]) & (x <= xnodes[2])
f2 = (xnodes[2] - x)/h[1]
#
Nfuncs.append( Piecewise( (f1, cond1), (f2, cond2), (0, True) ) )

# Third
cond1 = (x >= xnodes[Nnodes-2]) & (x <= xnodes[Nnodes-1])
f1 = (x - xnodes[Nnodes-2])/h[Nnodes-2]
Nfuncs.append( Piecewise( (f1, cond1), (0, True) ) )


Lnum = 2.0

# import numpy as np
# import matplotlib.pyplot as plt
# NptsPlot = 200
# xplt = np.linspace(0.0, Lnum, 200)
# yplt = np.zeros(NptsPlot)
# plt.clf()
# for ibasis in range(Nnodes):
#     for i in range(NptsPlot):
#         yplt[i] = Nfuncs[ibasis].subs({x: xplt[i], L: Lnum})
#     plt.plot(xplt, yplt, label="Basis-" + str(ibasis+1))
# plt.legend()
# plt.grid(True)
# plt.savefig("IMG_basis_lin_01.png", dpi=150)


# Functions, symbolic
T = Function("T")(x)
w = Function("w")(x)
k = Symbol("k", real=True, positive=True)
Q = Symbol("Q", real=True)



Tnodal = []
for i in range(Nnodes):
    Tnodal.append(Symbol("T_" + str(i), real=True))

T_expansion = 0
for i in range(Nnodes):
    T_expansion += Tnodal[i]*Nfuncs[i]

# Some concrete values
pprint(T_expansion.subs({x: 0.1, L: Lnum}))


# Stiffness term
expr1 = k * Derivative(w,x) * Derivative(T, x)
expr1s = expr1.subs({T: T_expansion, w: Nfuncs[0]})
pprint(expr1s)

expr1ss = expr1s.doit()
pprint(expr1ss)

I1 = Integral(expr1ss, (x,0,L))
I1

term1_first = I1.doit()


# Source term
expr2 = Integral( Q*w, (x,0,L) )
term2_first = expr2.subs({w: Nfuncs[0]}).doit()
term2_first


expr3 = k*w*Derivative(T,x)
expr3

term3_x0 = expr3.subs({ w: Nfuncs[0], x: 0, L: Lnum})
term3_x0


term3_xL = expr3.subs({ w: Nfuncs[0], x: Lnum, L: Lnum})
term3_xL

term3_first = term3_xL - term3_x0
term3_first

#

eq_first = term1_first - term2_first - term3_first
eq_first
q0 = Symbol("q_0", real=True)
term3_first
term3_first.args
term3_first.args[2]
#eq_first = eq_first.subs( {k*term3_first.args[2]: -q0} )
#eq_first

