using Printf
using LinearAlgebra
using SparseArrays

import PyPlot as plt

function calc_analytic_solution!( κ, Lx, Ly, t, xgrid, ygrid, Texact )
    Nterms = 100 # HARDCODED
    x = xgrid .- Lx/2
    y = ygrid .- Ly/2
    l = Lx/2
    b = Ly/2
    Nx = size(x,1) # XXX FIXME
    Ny = size(y,1)
    sumx = zeros(Ny, Nx)
    sumy = zeros(Ny, Nx)
    for m in 0:Nterms
        eterm = exp( -κ*(2*m+1)^2*pi^2*t/(4*l^2) )
        costerm = cos.( (2*m+1)*pi*x/(2*l) )
        sumx = sumx + (-1)^m/(2*m+1)*eterm*costerm
        eterm = exp(-κ*(2*m+1)^2*pi^2*t/(4*b^2))
        costerm = cos.((2*m+1)*pi*y/(2*b))
        sumy = sumy + (-1)^m/(2*m+1)*eterm*costerm
    end
    Texact[:,:] = (4/pi)^2 * sumx[:,:] .* sumy[:,:] 
    return
end


seconds_per_yr = 60*60*24*365 # seconds in 1 year
# physical parameters
κ  = 1.0e-6    # thermal diffusivity, m^2/s
Lx = 2.0e3   # width of domain, m
Ly = 1.0e3     # depth of domain, m
H  = 0*1e-9  # heat source, o K/s

Tb = 0.0 # fixed boundary temperature
Ti = 1.0 # initial temperature

Ndim = 2
NelsX = 100
NelsY = 100
Nelements = NelsX*NelsY
#
NnodesPerElement = 4
NnodesX = NelsX + 1
NnodesY = NelsY + 1
NnodesTotal = NnodesX*NnodesY

dx = Lx/NelsX   # element length in x-direction
dy = Ly/NelsY   # element length in y-direction
    
dt = 10*seconds_per_yr # time step (s)
Ntime = 1  # number of time steps to perform
    
Ktensor = diagm( 0 => κ*ones(2) ) # thermal diffusivity tensor

g_coord = zeros(Ndim,NnodesTotal)
ip = 0
# y first
for i in 1:NnodesX
    for j in 1:NnodesY
        ip = ip + 1
        g_coord[1,ip] = (i-1)*dx
        g_coord[2,ip] = (j-1)*dy
    end
end

# Connectivity matrix
g_num = zeros(Int64, NnodesPerElement,Nelements) 
# grid of global node numbers
gnumbers = reshape( collect(1:NnodesTotal), NnodesY, NnodesX )
iel = 0 # element counter
# y first
for i in 1:NelsX
    for j in 1:NelsY
        iel = iel + 1
        # clockwise numbering
        g_num[1,iel] = gnumbers[j,i] # local node 1
        g_num[2,iel] = gnumbers[j+1,i] # local node 2
        g_num[3,iel] = gnumbers[j+1,i+1] # local node 3
        g_num[4,iel] = gnumbers[j,i+1] # local node 4
    end
end

# find boundary nodes
SMALL_XY = 0.01*min(dx,dy)

# small fraction of dx or dy
bx0 = findall(g_coord[1,:] .<= 0.0 + SMALL_XY)  # nodes on x=0 boundary
bxn = findall(g_coord[1,:] .>= Lx  - SMALL_XY)  # nodes on x=Lx boundary
by0 = findall(g_coord[2,:] .<= 0.0 + SMALL_XY)  # nodes on y=0 boundary
byn = findall(g_coord[2,:] .>= Ly  - SMALL_XY)  # nodes on y=Ly boundary

# define boundary conditions
bcdof = unique(vcat(bx0, bxn, by0, byn))  # boundary nodes
Nbcdof = length(bcdof)
bcval = Tb*ones(Nbcdof)  # boundary temperatures
println("Nbcdof = ", Nbcdof)

NintegPoints = 4
# gauss integration data
integPoints = zeros(NintegPoints,Ndim) # location of points
root3 = 1.0/sqrt(3)
integPoints[1,1] = -root3; integPoints[1,2] =  root3;
integPoints[2,1] =  root3; integPoints[2,2] =  root3;
integPoints[3,1] = -root3; integPoints[3,2] = -root3;
integPoints[4,1] =  root3; integPoints[4,2] = -root3;
wIntegPoints = ones(4) # weights
# save shape functions and their derivatives in local coordinates
# evaluated at integration points
fun_s = zeros(NintegPoints,NintegPoints)
der_s = zeros(Ndim,NintegPoints,NintegPoints)
der = zeros(Ndim,NintegPoints)
#
for k in 1:NintegPoints
    #
    ξ = integPoints[k,1]
    η = integPoints[k,2]
    #
    ηm = 0.25*(1.0 - η)
    ηp = 0.25*(1.0 + η)
    ξm = 0.25*(1.0 - ξ)
    ξp = 0.25*(1.0 + ξ)
    #
    fun = 4*[ ξm*ηm, ξm*ηp, ξp*ηp, ξp*ηm ]
    fun_s[k,:] = fun  # shape functions
    der[1,1] = -ηm
    der[1,2] = -ηp
    der[1,3] =  ηp
    der[1,4] =  ηm
    #
    der[2,1] = -ξm
    der[2,2] =  ξm
    der[2,3] =  ξp
    der[2,4] = -ξp
    der_s[:,:,k] = der # derivatives of shape function
end

#for i in 1:NintegPoints
#    println("\ni = ", i)
#    display(der_s[:,:,i]); println()
#end


# matrix integration and assembly
KM = zeros(NnodesPerElement,NnodesPerElement)
MM = zeros(NnodesPerElement,NnodesPerElement)
F = zeros(NnodesPerElement)
# Matrices and vectors
ff = zeros(Float64, NnodesTotal)
LHS = spzeros(Float64, NnodesTotal, NnodesTotal)
RHS = spzeros(Float64, NnodesTotal, NnodesTotal)
for iel in 1:Nelements # sum over elements
    num = g_num[:,iel] # element nodes
    coord = g_coord[:,num]' # element coordinates
    KM[:] .= 0.0
    MM[:] .= 0.0
    F[:] .= 0.0
    #
    for k in 1:NintegPoints
        fun = fun_s[k,:]  # shape functions
        der = der_s[:,:,k] # der. of shape functions in local coordinates
        jac = der*coord # jacobian matrix
        detjac = det(jac) # det. of jacobian
        invjac = inv(jac) # inv. of jacobian
        deriv = invjac*der # der. of shape fun. in physical coords.
        #
        KM = KM + deriv'*Ktensor*deriv*detjac*wIntegPoints[k]
        MM = MM + fun*fun'*detjac*wIntegPoints[k]
        F = F + fun*H*detjac*wIntegPoints[k]
    end
    # assemble global matrices and vector
    LHS[num,num] = LHS[num,num] + MM/dt + KM
    RHS[num,num] = RHS[num,num] + MM/dt
    ff[num] = ff[num] + F
end

Tnum = Ti*ones(NnodesTotal) # initial condition
b = zeros(NnodesTotal)
t = 0.0
for itime in 1:Ntime
    t = t + dt # update time
    #println(displ[100])
    b[:] = RHS*Tnum + ff # form rhs global vector
    #
    LHS[bcdof,:] .= 0.0 # zero the boundary equations
    for i in bcdof
        LHS[i,i] = 1.0
    end
    b[bcdof] = bcval # set boundary values
    #
    # Solve linear system
    Tnum[:] = LHS \ b 
    #ldiv!(displ, factorLHS, b)
end

println("Last t = ", t)
xgrid = reshape(g_coord[1,:],NnodesY,NnodesX)
ygrid = reshape(g_coord[2,:],NnodesY,NnodesX)
Texact = zeros(NnodesY,NnodesX)
calc_analytic_solution!( κ, Lx, Ly, t, xgrid, ygrid, Texact )
Tnumr = reshape(Tnum, (NnodesY,NnodesX))
diffT = sum(abs.( Texact .- Tnumr ))/(NnodesTotal)
println("MAE = ", diffT)
