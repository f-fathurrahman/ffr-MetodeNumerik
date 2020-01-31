using Printf

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("heat_1d_euler_imp.jl")

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

function plot_all(u, x, t; prefix="IMG_")
    
    Nt = size(t,1) - 1
    
    println("Plotting the solution, please wait ...")    
    
    for i in 1:Nt+1
        plt.clf()
        label_t = @sprintf("t=%f", t[i])
        println(label_t)
        plt.plot(x, u[:,i], marker="o", label=label_t)
        plt.plot(x, sol_01_analitik.(x, t[i]), label="exact")
        plt.ylim(0.0, 1.1)
        plt.legend()
        filename = @sprintf("IMG_ex01_heat_1d_euler_imp_%04d.png", i)
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

    u_imp, x_imp, t_imp = heat_1d_euler_imp( α, xf, tf, initial_temp, bx0, bxf, Nx, Nt )

    #plt.clf()
    #lbl_str = @sprintf("t=%f", t_imp[1])
    #plt.plot(x_imp, u_imp[:,1], marker="o", label=lbl_str)
    #lbl_str = @sprintf("t=%f", t_imp[end])
    #plt.plot(x_imp, u_imp[:,end], marker="o", label=lbl_str)
    #plt.legend()
    #plt.savefig("IMG_ex01_heat_1d_euler_imp.png", dpi=150)

    plot_all( u_imp, x_imp, t_imp )

end

main()