using Printf
using SpecialFunctions: erf

import PyPlot
const plt = PyPlot

function main()

    # Parameters
    T_L = 10.0
    q = 0.0
    k = 1.0
    L = 1.0
    alpha = 0.1
    A = 10.0
    
    function Q(x)
        return A*exp( -alpha*(x-L/2)^2 )
    end

    function T_analytic(x)
        out1 = T_L + q*(L - x)/k + A*( (3*sqrt(pi)*L*erf(L*sqrt(alpha)/2)/(4*sqrt(alpha)) - sqrt(pi)*L*erf(L*sqrt(alpha)/2 - sqrt(alpha)*x)/(4*sqrt(alpha)) - exp(-L^2*alpha/4 + L*alpha*x - alpha*x^2)/(2*alpha) + exp(-L^2*alpha/4)/(2*alpha) - sqrt(pi)*x*erf(L*sqrt(alpha)/2)/(2*sqrt(alpha)) + sqrt(pi)*x*erf(L*sqrt(alpha)/2 - sqrt(alpha)*x)/(2*sqrt(alpha)))/k )
        return out1
    end

    function q_L()
        out1 = -sqrt(pi)*A.*erf(L.*sqrt(alpha)/2)./sqrt(alpha) - q
        return out1
    end

#=
    Nelements = 4
    Nnodes = Nelements + 1
    NnodesPerElement = 2 # linear element
=#
    Nnodes = 50
    x = range(0, stop=L, length=Nnodes) # nodes

#=
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
    # Boundary term
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
=#

    T_exact = T_analytic.(x)
    Q_plt = Q.(x)
    #display(T_exact); println()
    plt.clf()
    plt.plot(x, T_exact, label="exact", marker="o")
    #plt.plot(x, Q_plt, label="Q_plt", marker="o")
    plt.legend()
    plt.grid()
    plt.savefig("IMG_ex01_gaussian_source.pdf")

    println("q_0 = ", q)
    println("T_0 = ", T_exact[1])

    println("q_L = ", q_L())
    println("T_L = ", T_L)
end

main()
