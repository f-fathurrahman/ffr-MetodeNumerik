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



function main(; do_plot=false)
    
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

    # Boundary condition
    Tb = 0.0
    bcdof = [1, NnodesTotal] # boundary nodes
    bcval = [Tb, Tb] # boundary values

    # Connectivity and equation numbering
    g_num = zeros(Int64,NnodesPerElement,Nelements);
    g_num[1,:] = collect(1:NnodesTotal-1)
    g_num[2,:] = collect(2:NnodesTotal)

    # Matrices and vectors
    ff = zeros(Float64, NnodesTotal)  # system load vector
    LHS = spzeros(Float64, NnodesTotal, NnodesTotal) # system lhs matrix
    RHS = spzeros(Float64, NnodesTotal, NnodesTotal) # system rhs matrix

    # ---------------
    # Matrix assembly
    # ---------------

    for iel = 1:Nelements
    
        num = g_num[:,iel] # get equation number

        dx = abs( g_coord[num[2]] - g_coord[num[1]] ) # length of the element
    
        # Mass matrix
        MM = dx*[ 1/3 1/6;
                  1/6 1/3]

        # Stiffness
        KM = [ κ/dx -κ/dx;
              -κ/dx  κ/dx]

        # Load vector
        F = dx*H*[1/2; 1/2]

        LHS[num,num] = LHS[num,num] + MM/dt + KM

        RHS[num,num] = RHS[num,num] + MM/dt

        ff[num] = ff[num] + F
    end

    Texact = zeros(NnodesTotal)
    Tnum = zeros(NnodesTotal)
    for i in 1:NnodesTotal
        Tnum[i] = sin(π*g_coord[i])
    end
    b = zeros(NnodesTotal)

    LHS[1,2] = 0.0
    LHS[1,1] = 1.0
    LHS[NnodesTotal,NnodesTotal-1] = 0.0
    LHS[NnodesTotal,NnodesTotal] = 1.0
    factorLHS = lu(LHS)

    plt.clf()
    plt.plot(g_coord, Tnum, marker="o")
    plt.grid(true)
    plt.savefig("IMG_case2_initial.png", dpi=150)

    # Time loop
    t = 0.0
    for ni in 1:Ntime

        t = t + dt

        # Boundary condition, pulled out from the time loop
        #LHS[1,2] = 0.0
        #LHS[1,1] = 1.0
        #LHS[NnodesTotal,NnodesTotal-1] = 0.0
        #LHS[NnodesTotal,NnodesTotal] = 1.0

        #println("After setting BC, RHS = ")
        #display(RHS); println()

        @views b[:] = RHS*Tnum + ff
        b[bcdof] .= bcval

        ldiv!(Tnum, factorLHS, b)
        #Tnum[:] = LHS\b
    end

    println("t = ", t)
    calc_analytic_solution!(a, t, g_coord, Texact)
    for i in 1:NnodesTotal
        @printf("%3d %18.10f %18.10f\n", i, Tnum[i], Texact[i])
    end
    rmse = sqrt( sum((Texact - Tnum).^2)/NnodesTotal )
    println("RMS error = ", rmse)

    plt.clf()
    plt.plot(g_coord, Tnum, label="numeric")
    plt.plot(g_coord, Texact, label="exact")
    plt.grid(true)
    plt.legend()
    plt.savefig("IMG_case2_final.png", dpi=150)

end

@time main()
#@time main()