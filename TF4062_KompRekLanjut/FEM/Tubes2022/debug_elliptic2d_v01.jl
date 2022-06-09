import MAT
using Printf
using LinearAlgebra
using SparseArrays

meshfile = "mesh_cylinder.mat"
vars = MAT.matread(meshfile)
g_coord = vars["g_coord"]
g_num = convert(Matrix{Int64}, vars["g_num"])

NnodesTotal = size(g_coord,2) # total number of nodes
Nelements = size(g_num,2) # total number of elements
Lx = maximum(g_coord[1,:]) # x-dimension of domain
Ly = maximum(g_coord[2,:]) # y-dimension of domain

println("NnodesTotal = ", NnodesTotal)
println("Nelements = ", Nelements)

# physical parameters
U = 1.0 # imposed boundary velocity (=phi/dx)

# numerical parameters
NnodesPerElement = 3 # number of nodes in 1 element
Ndim = 2  # number of spatial dimensions in problem

# gauss integration data
NintegPoints = 3  # number of Gauss integration points
points = zeros(NintegPoints, Ndim)
points[1,1] = 0.5; points[1,2] = 0.5
points[2,1] = 0.5; points[2,2] = 0.0
points[3,1] = 0.0; points[3,2] = 0.5
c = 0.5  # triangle factor
wts = c*(1/3)*ones(NintegPoints) # weights

# save shape functions and their derivatives in local coordinates
# evaluated at integration points
fun_s = zeros(NintegPoints,NintegPoints)
der_s = zeros(Ndim,NintegPoints,NintegPoints)
der = zeros(Ndim,NintegPoints)    
for k in 1:NintegPoints
    L1 = points[k,1]
    L2 = points[k,2]
    L3 = 1.0 - L1 - L2
    fun = [L1, L2, L3]
    fun_s[k,:] = fun # shape functions
    der[1,1] = 1.0; der[1,2] = 0.0; der[1,3] = -1.0
    der[2,1] = 0.0; der[2,2] = 1.0; der[2,3] = -1.0
    der_s[:,:,k] = der # derivative of shape function
end


# data required to impose Neumann boundary conditions
# establish elements (with 2 nodes) on x=0 and x=Lx boundaries
# save the local node indices of the 2 boundary nodes
SMALL = 0.001
IDXTRUE = [false,false,false]
IDXNODES = [1,2,3]
Nx0found = 0
for iel in 1:Nelements
    ix1 = g_num[1,iel]
    ix2 = g_num[2,iel]
    ix3 = g_num[3,iel]
    IDXTRUE[1] = g_coord[1,ix1] <= (0.0 + SMALL)
    IDXTRUE[2] = g_coord[1,ix2] <= (0.0 + SMALL)
    IDXTRUE[3] = g_coord[1,ix3] <= (0.0 + SMALL)
    if count(IDXTRUE) == 2
        Nx0found = Nx0found + 1
    end
end

belx0 = zeros(Int64,Nx0found)
iix0 = zeros(Int64,2,Nx0found)
ii0 = 0
for iel in 1:Nelements
    ix1 = g_num[1,iel]
    ix2 = g_num[2,iel]
    ix3 = g_num[3,iel]
    IDXTRUE[1] = g_coord[1,ix1] <= (0.0 + SMALL)
    IDXTRUE[2] = g_coord[1,ix2] <= (0.0 + SMALL)
    IDXTRUE[3] = g_coord[1,ix3] <= (0.0 + SMALL)
    if count(IDXTRUE) == 2 # 2 nodes on x=0
        ii0 = ii0 + 1
        belx0[ii0] = iel
        iix0[:,ii0] .= IDXNODES[IDXTRUE]
    end
end
println(iix0)

Nxnfound = 0
for iel in 1:Nelements
    ix1 = g_num[1,iel]
    ix2 = g_num[2,iel]
    ix3 = g_num[3,iel]
    IDXTRUE[1] = g_coord[1,ix1] >= (Lx - SMALL)
    IDXTRUE[2] = g_coord[1,ix2] >= (Lx - SMALL)
    IDXTRUE[3] = g_coord[1,ix3] >= (Lx - SMALL)
    if count(IDXTRUE) == 2
        Nxnfound = Nxnfound + 1
    end
end

belxn = zeros(Int64,Nxnfound)
iixn = zeros(Int64,2,Nxnfound)
iin = 0
for iel in 1:Nelements
    ix1 = g_num[1,iel]
    ix2 = g_num[2,iel]
    ix3 = g_num[3,iel]
    IDXTRUE[1] = g_coord[1,ix1] >= (Lx - SMALL)
    IDXTRUE[2] = g_coord[1,ix2] >= (Lx - SMALL)
    IDXTRUE[3] = g_coord[1,ix3] >= (Lx - SMALL)
    if count(IDXTRUE) == 2
        iin = iin + 1
        belxn[iin] = iel
        iixn[:,iin] .= IDXNODES[IDXTRUE]
    end
end
#println(iixn)

#println("Nx0found = ", Nx0found)
#println("Nxnfound = ", Nxnfound)

# global matrices and vectors
bv = zeros(NnodesTotal) # global rhs vector
LHS = spzeros(NnodesTotal,NnodesTotal) # global lhs matrix
    
#----------------------------------------------------
# matrix integration and assembly
#----------------------------------------------------
KM = zeros(NnodesPerElement,NnodesPerElement) # initialise element stiffness matrix
for iel in 1:Nelements # loop over elements
    num = .... # element nodes
    coord = .... # element coordinates
    KM .= 0.0
    for k in 1:NintegPoints # integration loop
        der = .... # der. of shape functions in local coordinates
        jac = ....
        # jacobian matrix
        detjac = .... # det. of jacobian
        invjac = .... # inv. of jacobian
        deriv = .... # der. of shape fun. in physical coords.
        KM += .... # stiffness matrix
    end
    # assemble global matrix
    for j in 1:NnodesPerElement, i in 1:NnodesPerElement
        ii = num[i]
        jj = num[j]
        LHS[ii,jj] += KM[i,j]
    end
end

# compute boundary vector for Neumann condition on x=0
mmv = zeros(Int64,NnodesPerElement)
for i in 1:size(belx0,1)
    iel = belx0[i]
    num = g_num[:,iel]
    y = g_coord[2,num]'
    ln = iix0[:,i]
    mmv .= 0
    mmv[ln] .= 1
    b = abs.( diff(y[ln]) )
    Fb = 0.5 * U * b .* mmv
    bv[num] = bv[num] - Fb
end

# compute boundary vector for Neumann condition on x=Lx
for i in 1:size(belxn,1)
    iel = belxn[i]
    num = g_num[:,iel]
    y = g_coord[2,num]
    ln = iixn[:,i]
    mmv .= 0
    mmv[ln] .= 1
    b = abs.( diff(y[ln]) )
    Fb = 0.5 * U * b .* mmv
    bv[num] = bv[num] + Fb
end

# NOTE: The solution is unique up to a constant.
# The result may be different from MATLAB or Octave
displ = ....


#
# Post-processing
#
F = zeros(NnodesPerElement) # Preallocate
U = zeros(NnodesTotal)
V = zeros(NnodesTotal)
ugp = zeros(NnodesPerElement,Nelements)
vgp = zeros(NnodesPerElement,Nelements)
for iel in 1:Nelements # loop over elements
    num = g_num[:,iel] # element nodes
    coord = g_coord[:,num]' # element coordinates
    fill!(KM, 0.0)
    fill!(F, 0.0) 
    for k in 1:NintegPoints
        fun = fun_s[k,:] # shape functions
        der = der_s[:,:,k] # der. of shape functions in local coordinates
        jac = der*coord   # jacobian matrix
        detjac = det(jac) # det. of jacobian
        invjac = inv(jac) # inv. of jacobian
        deriv = invjac*der # der. of shape fun. in physical coords.
        #
        # Calculate u and v
        u = dot(deriv[1,:], displ[num])
        v = dot(deriv[2,:], displ[num])
        ugp[k,iel] = u # save integration point velocity
        vgp[k,iel] = v # save integration point velocity
    end
end



#
# Export data ke VTK (dapat divisualisasikan dengan Paraview)
#
using WriteVTK: MeshCell, VTKCellTypes, vtk_grid, VTKPointData, VTKCellData
points = g_coord
Ncells = size(g_num,2)
cells = Vector{MeshCell}(undef,Ncells)
for i in 1:Ncells
    @views cells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, g_num[:,i])
end

cell_velocity = zeros(2,Nelements)
cell_velocity[1,:] = ugp[1,:] # ambil nilai pada node (lokal) pertama
cell_velocity[2,:] = vgp[1,:] # ambil nilai pada node (lokal) pertama
vtk_grid("TEMP_elliptic", points, cells) do vtk
    vtk["potential", VTKPointData()] = displ
    vtk["cell_velocity", VTKCellData()] = cell_velocity
end
