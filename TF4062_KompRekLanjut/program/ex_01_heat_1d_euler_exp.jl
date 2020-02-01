using Printf

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("heat_1d_euler_exp.jl")

function initial_temp(x)
    return sin(pi*x)
end

# Syarat batas kiri
function bx0( t )
    return 0.0
end

# Syarat batas kanan
function bxf( t )
    return 0.0
end

# Solusi analitik
function sol_01_analitik(x, t)
    return sin(pi*x) * exp(-pi^2 * t)
end

function plot_all(u, x, t)

    Nt = size(t,1)
    
    println("Plotting the solution, please wait ...")    
    
    for i in 1:Nt
        plt.clf()
        label_t = @sprintf("t=%f", t[i])
        println(label_t)
        plt.plot(x, u[:,i], marker="o", label=label_t)
        plt.plot(x, sol_01_analitik.(x, t[i]), label="exact")
        plt.ylim(0.0, 1.1)
        plt.legend()
        filename = @sprintf("IMG_ex01_heat_1d_euler_exp_%04d.png", i)
        plt.savefig(filename, dpi=150)
    end
    
    println("Done")

end


function main()
    α = 1.0
    xf = 1.0
    tf = 0.1
    Nx = 25
    Nt = 200

    u_exp, x_exp, t_exp = heat_1d_euler_exp( α, xf, tf, initial_temp, bx0, bxf, Nx, Nt )

    plot_all( u_exp, x_exp, t_exp )
end

main()