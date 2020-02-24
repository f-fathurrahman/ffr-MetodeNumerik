using Printf

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

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

function plot_all(u, x, t; prefix="IMG_")
    
    Nt = size(t,1)
    
    println("Plotting the solution, please wait ...")    
    
    for i in 1:Nt
        plt.clf()
        label_t = @sprintf("t=%f", t[i])
        println(label_t)
        plt.plot(x, u[:,i], marker="o", label=label_t)
        plt.plot(x, analytic_solution.(x, t[i]), label="exact")
        plt.ylim(0.0, 1.1)
        plt.legend()
        filename = @sprintf("IMG_ex01_diffusion_1d_CN_%04d.png", i)
        plt.savefig(filename, dpi=150)
    end
    
    println("Done")

end


function main()
    α = 1.0
    L = 1.0
    T = 0.2
    Nx = 25
    Nt = 400

    u, x, t = diffusion_1d_CN( L, Nx, T, Nt, α, initial_temp, bx0, bxf, source_term )

    #plot_all( u, x, t )
    u_a = analytic_solution.(x, t[end])
    u_n = u[:,end]
    rmse = sqrt( sum((u_a - u_n).^2)/Nx )
    mean_abs_diff = sum( abs.(u_a - u_n) )/Nx
    @printf("RMS error            = %15.10e\n", rmse)
    @printf("Means abs diff error = %15.10e\n", mean_abs_diff)

end

main()