import PyPlot
const plt = PyPlot

function initial_cond(x; A=1.0, L=1.0, σ=0.02)
    return A*exp( -0.5*( (x - L/10)/σ )^2 )
end


function initial_cond_v2(x; A=1.0, L=1.0, σ=0.02)
    if x < L/5
        return A*cos(5*pi/L*(x-L/10))
    else
        return 0.0
    end
end

L = 1.0
v = 1.0
A = 1.0
σ = 0.02

Δx = 0.04
x = collect(range(0.0, stop=1.0, step=Δx))
#u = initial_cond.(x; A=A, L=L, σ=σ)
u = initial_cond_v2.(x; A=A, L=L, σ=σ)

Δt = 0.001
C = v*Δt/Δx

Nx = length(x)
unp1 = zeros(Float64, Nx)

plt.clf()
plt.plot(x, u, label="u")
plt.ylim(-0.1, 1.1)
plt.grid(true)
plt.legend()
plt.savefig("IMG_adv_" * string(0) * ".png", dpi=150)

INTERVAL_DO_PLOT = 50

for i in 1:1000
    println("Time = ", i*Δt)
    unp1[1] = 0.0 # Boundary condition
    unp1[Nx] = 0.0 # Boundary condition
    for i in 2:Nx-1
        unp1[i] = u[i] - 0.5*C*( u[i+1] - u[i-1] )
    end
    @views u[:] = unp1[:] # copy new value to old value
    #
    if (i % INTERVAL_DO_PLOT == 0) || (i == 1)
        plt.clf()
        plt.plot(x, u, label="u")
        plt.ylim(-0.1, 1.1)
        plt.grid(true)
        plt.legend()
        plt.savefig("IMG_adv_" * string(i) * ".png", dpi=150)
    end
end

