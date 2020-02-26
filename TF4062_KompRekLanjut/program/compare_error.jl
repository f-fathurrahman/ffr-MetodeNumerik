using Printf

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("diffusion_1d_explicit.jl")
include("diffusion_1d_implicit.jl")
include("diffusion_1d_CN.jl")

function initial_temp(x)
    return sin(π*x)
end

function bx0( t )
    return 0.0
end

function bxf( t )
    return 0.0
end

function source_term(x, t)
    return 0.0
end

function analytic_solution(x, t)
    return sin(π*x) * exp(-π^2 * t)
end


function calc_error(Nx, Nt, solver)
    
    α = 1.0
    L = 1.0
    T = 0.2

    Δx = L/(Nx-1)
    Δt = T/(Nt-1)

    u, x, t = solver( L, Nx, T, Nt, α, initial_temp, bx0, bxf, source_term )

    u_a = analytic_solution.(x, t[end])
    u_n = u[:,end]
    rmse = sqrt( sum((u_a - u_n).^2)/Nx )
    mae = sum( abs.(u_a - u_n) )/Nx

    return Δx, Δt, rmse, mae

end

function main()
    Nt = 500
    for Nx in 10:10:100
        #Δx, Δt, rmse, mae = calc_error( Nx, Nt, diffusion_1d_implicit )
        #Δx, Δt, rmse, mae = calc_error( Nx, Nt, diffusion_1d_CN )
        Δx, Δt, rmse, mae = calc_error( Nx, Nt, diffusion_1d_explicit )
        @printf("%18.10f %18.10e %18.10e\n", Δx, rmse, mae)
    end
end

main()