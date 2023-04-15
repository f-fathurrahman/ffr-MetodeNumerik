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
    f1 = (1 - 3*ξ) * (1 - 3*ξ/2) * (1 - ξ)
    f2 = 9*ξ * (1 - 3*ξ/2) * (1 - ξ)
    f3 = -9*ξ/2 * (1 - 3*ξ) * (1 - ξ)
    f4 = ξ * (1 - 3*ξ) * (1 - 3*ξ/2)
    return f1, f2, f3, f4

# Spatial variable
x = Symbol("x", real=True)

# Domain length
#L = Symbol("L", real=True, positive=True)
L = 1.0

# Setup cubic element
Nelements = 2
NnodesPerElement = 4
Nnodes = (NnodesPerElement-1)*Nelements + 1
#
# 1 2 3 4 : 1 element -> 4 (3 + 1)
# 4 5 6 7 : 2 element -> 7 (2*3 + 1)
# 7 8 9 10: 3 element -> 10 (3*3 + 1)
# 10 11 12 13: 4 element -> 13 (4*3 + 1)



# Mapping between element index to (global) nodes index
elm2nodes = []
elm2nodes.append([0,1,2,3])
for iel in range(1,Nelements):
    i0 = elm2nodes[iel-1][0] + 3
    i1 = elm2nodes[iel-1][1] + 3
    i2 = elm2nodes[iel-1][2] + 3
    i3 = elm2nodes[iel-1][3] + 3
    elm2nodes.append([i0,i1,i2,i3])


# This loop is not efficient
nodes2elm = []
for inode in range(Nnodes):
    idx_node = []
    for iel in range(Nelements):
        il, im1, im2, ir = elm2nodes[iel]
        if inode == il:
            idx_node.append(iel)
        if inode == im1:
            idx_node.append(iel)
        if inode == im2:
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
for iel in range(Nelements):
    il = elm2nodes[iel][0]
    ir = elm2nodes[iel][3]
    h.append(xnodes[ir]-xnodes[il])


assert(len(xnodes) == Nnodes)
assert(len(h) == Nelements)

funcs_elm = []
conds_elm = [] 
for iel in range(Nelements):
    il, im1, im2, ir = elm2nodes[iel]
    xl = xnodes[il]
    xm1 = xnodes[im1]
    xm2 = xnodes[im2]
    xr = xnodes[ir]
    ξ, cond_elm = my_mapping(x, xl, xr, L)
    f1, f2, f3, f4 = gen_local_basis(ξ)
    funcs_elm.append([f1, f2, f3, f4])
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
        idx_local_node = idx_local_node % 4
        print("idx_local_node = ", idx_local_node)
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
        # In 1d we only have at max 2 common elements
        raise RuntimeException("Not supported")

