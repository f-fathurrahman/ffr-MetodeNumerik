using Printf

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("diffusion_1d_explicit.jl")

const L =  1.0
const α = 1.0

function analytic_solution(x, t)
    return 5*t*x*(L - x)
end

function source_term(x, t)
    return 10*α*t + 5*x*(L - x)
end

function initial_cond(x)
    return analytic_solution(x, 0.0)
end

function bx0(t)
    return 0.0
end

function bxL(t)
    return 0.0
end

function main()
    T = 0.1
    Nx = 21
    F = 0.5
    dx = L/(Nx-1)
    Δt = F*dx^2/α
    println(T/Δt)
    Nt = round(Int64,T/Δt) + 1
    println("Nt = ", Nt)

    u, x, t = diffusion_1d_explicit(L, Nx, T, Nt, α, initial_cond, bx0, bxL, source_term)

    u_e = analytic_solution.(x, t[end])
    diff_u = maximum(abs.(u_e - u[:,end]))
    println("diff_u = ", diff_u)

    #for n = 1:Nt
    #    plt.clf()
    #    plt.plot(x, u[:,n], marker="o")
    #    plt.ylim(0.0, 0.13)
    #    filename = @sprintf("IMG_ver_%04d.png", n)
    #    println(filename)
    #    plt.savefig(filename, dpi=150)
    #end

end

main()