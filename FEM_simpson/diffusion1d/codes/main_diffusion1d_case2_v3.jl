# Based on: Simpson - Practical Finite Element Modelling in Earth Science

using Printf
using LinearAlgebra
using SparseArrays
import PyPlot as plt

#=
a^2 ∂^2/∂x^2 u(x,t) = ∂/∂t u(x,t)
u(x,0) = sin(πx) 
u(0,t) = 0
u(1,t) = 0

Analytic solution: u(x,t) = exp(-a^2 * π^2 * t) * sin(πx)
=#


function calc_analytic_solution!(
    a::Float64, t::Float64, x,
    Texact::Vector{Float64}
)
    fill!(Texact, 0.0)
    Npoints = size(x,1)
    for ip in 1:Npoints
        Texact[ip] = exp(-a^2 * π^2 * t) * sin(π*x[ip])
    end
    return
end

    
# -------------------
# Physical parameters
# -------------------
a2 = 2.0
κ = a2
a = sqrt(a2)
H = 0.0  # Heat source term (K/s)

# Spatial
Lx = 1.0

# Time
dt = 0.0001
Ntime = 2000

# Basis
Nelements = 300  # Total number of elements
NnodesPerElement = 2  # no. of nodes per element
NnodesTotal = Nelements + 1  # Total number of nodes
dx = Lx/Nelements

# global grid coordinates
g_coord = 0.0:dx:Lx  # the size NnodesTotal

# Integration stuffs
NIntegPoints = 2
integPoints = [-sqrt(1/3), sqrt(1/3)]
wIntegPoints = [1.0, 1.0]

# Shape functions and their derivatives in local coordinates
# 2nd dimension = number of nodes per element
fun_s = zeros(NIntegPoints,2)
der_s = zeros(NIntegPoints,2)
for k in 1:NIntegPoints
    fun_s[k,:] = [ (1-integPoints[k])/2, (1+integPoints[k])/2 ]
    der_s[k,:] = [-0.5, 0.5]
end

# Boundary condition
Tb = 0.0
bcdof = [1, NnodesTotal] # boundary nodes
bcval = [Tb, Tb] # boundary values

# Connectivity and equation numbering
g_num = zeros(Int64,NnodesPerElement,Nelements);
g_num[1,:] = collect(1:NnodesTotal-1)
g_num[2,:] = collect(2:NnodesTotal)


# Local matrices
MM = zeros(NnodesPerElement,NnodesPerElement)
KM = zeros(NnodesPerElement,NnodesPerElement)
F = zeros(NnodesPerElement)
# Global matrices and vectors
ff = zeros(Float64, NnodesTotal)  # system load vector
LHS = spzeros(Float64, NnodesTotal, NnodesTotal) # system lhs matrix
RHS = spzeros(Float64, NnodesTotal, NnodesTotal) # system rhs matrix
#
# Begin matrix assembly, loop over elements
#
for iel in 1:Nelements
    
    num = g_num[:,iel] # get equation number

    # length of the element
    dx = abs( g_coord[num[2]] - g_coord[num[1]] )
    
    KM[:] .= 0.0
    MM[:] .= 0.0
    F[:] .= 0.0

    for k in 1:NIntegPoints
        fun = fun_s[k,:]
        der = der_s[k,:]
        detjac = dx/2
        invjac = 2/dx
        deriv = der*invjac
        MM = MM + fun*fun' * detjac * wIntegPoints[k]
        KM = KM + deriv*deriv' * κ * detjac * wIntegPoints[k]
        F = F + H * fun * detjac * wIntegPoints[k]
    end

    LHS[num,num] = LHS[num,num] + MM/dt + KM
    RHS[num,num] = RHS[num,num] + MM/dt
    ff[num] = ff[num] + F
end

# Initial condition
Tnum = zeros(NnodesTotal)
for i in 1:NnodesTotal
    Tnum[i] = sin(π*g_coord[i])
end
b = zeros(NnodesTotal)

# Time loop
t = 0.0
for tstep in 1:Ntime
    t = t + dt # increment time
    @views b[:] = RHS*Tnum + ff
    # Set up Dirichlet BC
    # The LHS matrix
    LHS[1,2] = 0.0
    LHS[1,1] = 1.0
    LHS[NnodesTotal,NnodesTotal-1] = 0.0
    LHS[NnodesTotal,NnodesTotal] = 1.0
    # b vector is set to T values given by BC
    b[bcdof] .= bcval
    # Solve the linear system
    Tnum[:] = LHS\b
end
println("t = ", t)

#
# Compare with analytical solutions
#
Texact = zeros(NnodesTotal)
calc_analytic_solution!(a, t, g_coord, Texact)
#for i in 1:NnodesTotal
#    @printf("%3d %18.10f %18.10f\n", i, Tnum[i], Texact[i])
#end
rmse = sqrt( sum((Texact - Tnum).^2)/NnodesTotal )
println("RMS error = ", rmse)

