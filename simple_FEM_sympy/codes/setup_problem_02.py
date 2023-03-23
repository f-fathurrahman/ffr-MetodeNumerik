# Linear element, arbitrary number of nodes

from sympy import *

init_printing(use_unicode=True)

# Spatial variable
x = Symbol("x", real=True)

# Domain length
L = Symbol("L", real=True, positive=True)


# Setup linear element
Nelements = 3
Nnodes = Nelements + 1

# Setup the grid (same spacing)
xnodes = []
h = []
Δx = L/(Nnodes-1)
for i in range(Nnodes):
    xnodes.append(Integer(0) + i*Δx)
for i in range(Nelements):
    h.append(xnodes[i+1]-xnodes[i])

assert(len(xnodes) == Nnodes)
assert(len(h) == Nelements)

# Setup basis functions
Nfuncs = []

# First
cond1 = (x >= xnodes[0]) & (x <= xnodes[1])
f1 = (xnodes[1] - x)/h[0]
Nfuncs.append( Piecewise( (f1, cond1), (0, True) ) )

for i in range(1,Nnodes-1):
    cond1 = (x >= xnodes[i-1]) & (x <= xnodes[i])
    f1 = ( x - xnodes[i-1])/h[i-1]
    #
    cond2 = (x >= xnodes[i]) & (x <= xnodes[i+1])
    f2 = (xnodes[i+1] - x)/h[i]
    #
    Nfuncs.append( Piecewise( (f1, cond1), (f2, cond2), (0, True) ) )

# Third
cond1 = (x >= xnodes[Nnodes-2]) & (x <= xnodes[Nnodes-1])
f1 = (x - xnodes[Nnodes-2])/h[Nnodes-2]
Nfuncs.append( Piecewise( (f1, cond1), (0, True) ) )

