#=
Solves:

-k * d2T/dx2 = Q
0 < x L

-k dT/dx = q   for x = 0
T = T_{L} for x = L 
=#

using Printf

import PyPlot
const plt = PyPlot

function main()

    # Parameters
    T_L = 10.0
    q = 5.0
    k = 4.0
    L = 1.0
    Q = 10.0

    T_analytic(x) = T_L + q/k * (L - x) + Q/(2*k) * (L^2 - x^2)

    Nelements = 30
    Nnodes = Nelements + 1
    NnodesPerElement = 2 # linear element

    x = range(0, stop=L, length=Nnodes) # nodes

    h = zeros(Float64,Nelements)
    for i in 1:Nelements
        h[i] = x[i+1] - x[i]
    end

    # Global stiffness matrix.
    # For simplicity, we are not using sparse matrix
    K = zeros(Float64,Nnodes,Nnodes)

    # Global load vector
    f = zeros(Float64,Nnodes)
    f_local = zeros(Float64,NnodesPerElement)

    # Local stiffness matrix, allocate the memory here
    K_local = zeros(Float64,NnodesPerElement,NnodesPerElement)

    # connectivity
    NodesForElement = zeros(Int64,NnodesPerElement,Nelements)
    for iel in 1:Nelements
        # Special case
        NodesForElement[1,iel] = iel
        NodesForElement[2,iel] = iel + 1
    end

    for iel in 1:Nelements
        K_local[1,1] = k/h[iel]
        K_local[1,2] = -k/h[iel]
        K_local[2,1] = -k/h[iel]
        K_local[2,2] = k/h[iel]
        idx = NodesForElement[:,iel]
        #
        K[idx,idx] = K[idx,idx] + K_local
        #
        f_local[1] = 0.5*Q*h[iel]
        f_local[2] = 0.5*Q*h[iel]
        f[idx] = f[idx] + f_local
    end
    # Apply boundary condition Boundary term
    f[1] = f[1] + q
    #
    f[Nnodes-1] = f[Nnodes-1] + T_L*k/h[Nelements]
    
    #display(K); println()
    #display(f); println()

    # Solution
    T = zeros(Float64,Nnodes)
    T[Nnodes] = T_L
    idx_solve = 1:Nnodes-1
    T[idx_solve] = K[idx_solve,idx_solve]\f[idx_solve]
    display(T); println()

    T_exact = T_analytic.(x)
    display(T_exact); println()
    #plt.clf()
    #plt.plot(x, T_exact, label="exact", marker="o")
    #plt.legend()
    #plt.grid()
    #plt.savefig("IMG_ex01.pdf")

    # Flux at x=L (from the last equation that is dropped)
    q_L = k/h[Nelements] * (T_L - T[Nnodes-1]) + Q/2*h[Nelements]

    println("q_0 = ", q)
    println("q_L = ", q_L)


    println("q_L analytic: ", -L*Q - q)
end

main()
