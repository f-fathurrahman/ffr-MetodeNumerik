# environment: GRIDAP

import GLMakie
GLMakie.set_theme!(GLMakie.theme_dark())

using Infiltrator

function do_plot(x, y, u)
    fig = GLMakie.surface(x, y, u)
    return fig
end

function main_burgers()
    Lx = 2.0
    Ly = 4.0
    Nx = 41
    Ny = 81
    Nt = 12

    c = 1
    Δx = Lx / (Nx - 1)
    Δy = Ly / (Ny - 1)
    ν = 0.01
    σ = 0.2
    Δt = σ * Δx

    x = range(0.0, stop=Lx, length=Nx)
    y = range(0.0, stop=Ly, length=Ny)

    u = zeros(Float64, Nx, Ny)
    v = zeros(Float64, Nx, Ny)
    un = zeros(Float64, Nx, Ny)
    vn = zeros(Float64, Nx, Ny)

    # initial conditions
    range_ix = Integer(0.5/Δx):Integer(1/Δx+1)
    range_iy = Integer(0.5/Δy):Integer(1/Δy+1)
    u[range_ix, range_iy] .= 2.0
    v[range_ix, range_iy] .= 2.0

    println("Pass here 1")

    @exfiltrate

end
