using Revise
using Infiltrator

includet("common.jl")
includet("region.jl")


function test_main()

    Lx = 2.0
    Ly = 1.0
    Nx = 5
    Ny = 5
    region = Region(Lx, Ly, Nx, Ny)

    # primary grid
    xc = (region.idx_ic .- 1) .* region.dx
    yc = (region.idx_jc .- 1) .* region.dy

    # midcell grid
    xm = (region.idx_ic .- 0.5) .* region.dx
    ym = (region.idx_jc .- 0.5) .* region.dy
    
    xx, yy = meshgrid(xm, ym)
    xx = xx'
    yy = yy'

    Nxm = region.Nxm
    Nym = region.Nym
    # initial velocity field (icas=1)
    u = zeros(Float64, Nxm, Nym)
    sca = zeros(Float64, Nxm, Nym)

    NSE_F_init_KH!(
        region,
        Lx, Ly, xc, ym,
        1, Ly/4, 20, 0.25, 0.5*Lx,
        u
    )
    NSE_F_init_KH!(
        region,
        Lx, Ly, xm, ym,
        1, Ly/4, 20, 0.00, 0.5*Lx,
        sca
    )

    @infiltrate

    #=
    Rey = 1000
    Pec = 1000
    Tstepmax = 210
    nprint = 10
    niso = 10
    =#
end

#test_main()