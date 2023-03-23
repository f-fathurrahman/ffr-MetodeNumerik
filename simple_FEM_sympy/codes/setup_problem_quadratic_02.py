# Linear element, arbitrary number of nodes

from sympy import *

init_printing(use_unicode=True)

def my_mapping(x, x1, x2, L):
    assert(x2 > x1)
    cond1 = x >= x1
    cond2 = x <= x2
    ξ = (x - x1) * (L/(x2 - x1))
    return ξ, cond1 & cond2

def gen_local_basis(ξ):
    f1 = 1 - 3*ξ + 2*ξ**2
    f2 = 4*ξ*(1 - ξ)
    f3 = ξ*(2*(ξ) - 1)
    return f1, f2, f3

# Spatial variable
x = Symbol("x", real=True)

# Domain length
#L = Symbol("L", real=True, positive=True)
L = 1.0

# Setup quadratic element
Nelements = 3
Nnodes = 2*Nelements + 1
NnodesPerElement = 3

# Mapping between element index to (global) nodes index
elm2nodes = []
elm2nodes.append([0,1,2])
for iel in range(1,Nelements):
    i0 = elm2nodes[iel-1][0] + 2
    i1 = elm2nodes[iel-1][1] + 2
    i2 = elm2nodes[iel-1][2] + 2
    elm2nodes.append([i0,i1,i2])


# This loop is not efficient
nodes2elm = []
for inode in range(Nnodes):
    idx_node = []
    for iel in range(Nelements):
        il, im, ir = elm2nodes[iel]
        if inode == il:
            idx_node.append(iel)
        if inode == im:
            idx_node.append(iel)
        if inode == ir:
            idx_node.append(iel)
    nodes2elm.append(idx_node)



# Setup the grid (same spacing)
Δx = L/(Nnodes-1)
xnodes = []

for i in range(Nnodes):
    xnodes.append(Integer(0) + i*Δx)

h = []
for i in range(Nelements):
    h.append(xnodes[i+1]-xnodes[i])

assert(len(xnodes) == Nnodes)
assert(len(h) == Nelements)

funcs_elm = []
conds_elm = [] 
for iel in range(Nelements):
    il, im, ir = elm2nodes[iel]
    xl = xnodes[il]
    xm = xnodes[im]
    xr = xnodes[ir]
    ξ, cond_elm = my_mapping(x, xl, xr, L)
    f1, f2, f3 = gen_local_basis(ξ)
    funcs_elm.append([f1, f2, f3])
    conds_elm.append(cond_elm)

idx_local_node = 0
funcs_nodes = []
conds_nodes = []
for inode in range(Nnodes):
    iels = nodes2elm[inode]
    #
    func_node = []
    cond_node = []
    for iel in iels:
        idx_local_node = idx_local_node % 3
        #print("idx_local_node = ", idx_local_node)
        func_node.append(funcs_elm[iel][idx_local_node])
        cond_node.append(conds_elm[iel])
        idx_local_node += 1
    funcs_nodes.append(func_node)
    conds_nodes.append(cond_node)

# Now we are ready to build global interpolation functions as 
# list of Piecewise function
Nfuncs = []
for inode in range(Nnodes):
    assert len(conds_nodes[inode]) == len(funcs_nodes[inode])
    if len(conds_nodes[inode]) == 1:
        c = conds_nodes[inode][0]
        f = funcs_nodes[inode][0] 
        f = Piecewise( (f, c), (0, True) )
        Nfuncs.append(f)
    elif len(conds_nodes[inode]) == 2:
        #
        c0 = conds_nodes[inode][0]
        f0 = funcs_nodes[inode][0]
        #
        c1 = conds_nodes[inode][1]
        f1 = funcs_nodes[inode][1] 
        f = Piecewise( (f0, c0), (f1, c1), (0, True) )
        Nfuncs.append(f)
    else:
        raise RuntimeException("Not supported")

