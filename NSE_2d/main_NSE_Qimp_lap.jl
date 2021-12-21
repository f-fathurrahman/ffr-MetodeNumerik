using Printf

include("common.jl")
include("NSE_calc_lap.jl")
include("NSE_norm_L2.jl")
include("NSE_ADI_init.jl")
include("NSE_ADI_step.jl")

function NSE_fsource(Lx, Ly, x , y)
     aa = 2π/Lx
     bb = 2π/Ly
     fs = (aa*aa + bb*bb) * ( cos.(aa*x) .* sin.(bb*y) )
     return fs
end

function NSE_fexact(Lx, Ly, x, y)
    aa = 2π/Lx
    bb = 2π/Ly
    fex = cos.(aa*x) .* sin.(bb*y)
    return fex
end

function main()

    # Input parameters
    Lx = 1.0
    Ly = 1.0
    Nx = 51
    Ny = 51

    # 2D grid variables

    Nxm = Nx - 1
    Nym = Ny - 1
    dx = Lx/Nxm
    dy = Ly/Nym

    ic = collect(1:Nxm)
    jc = collect(1:Nym)

    xc = (ic .- 1)*dx
    xm = (ic .- 0.5)*dx
    
    yc = (jc .- 1)*dy
    ym = (jc .- 0.5)*dy

    ip = ic .+ 1
    ip[Nxm] = 1
    
    jp = jc .+ 1
    jp[Nym] = 1
    
    imm = ic .- 1
    imm[1] = Nxm
    
    jm = jc .- 1
    jm[1] = Nym

    xx, yy = meshgrid(xc, ym)
    xx = xx'
    yy = yy'

    #  Initialization
    u = zeros(Float64, Nxm, Nym)
    du = zeros(Float64, Nxm, Nym)
    hc = zeros(Float64, Nxm, Nym)
    rhs = zeros(Float64, Nxm, Nym)

    # Time step                   
    dt = 0.5/(1/dx^2 + 1/dy^2 )*100
    bx = 0.5*dt/dx^2
    by = 0.5*dt/dy^2

    amix, apix, alphx, xs2x = NSE_ADI_init(-bx*ones(Nxm), (1+2*bx)*ones(Nxm), -bx*ones(Nxm))
    amiy, apiy, alphy, xs2y = NSE_ADI_init(-by*ones(Nym), (1+2*by)*ones(Nym), -by*ones(Nym))
    
    # Time loop
    CONV = 1
    NiterMax = 10_000
    Niter = 0
    temps = 0

    while( (CONV > 1e-6) && (Niter <= NiterMax) )
        
        Niter = Niter + 1
        temps = temps + dt
        
        # RHS term
        @views rhs[:,:] = -0.5*dt*hc[:,:]
        @views hc[:,:]  = NSE_fsource(Lx, Ly, xx, yy)
        @views rhs[:,:] = rhs[:,:] + 1.5*dt*hc[:,:] + dt*NSE_calc_lap(u, dx, dy, imm, ip, jp, jm, ic, jc)
        
        # first step of ADI
        rhs = NSE_ADI_step(amix, apix, alphx, xs2x, rhs)
        # second step of ADI
        rhs = NSE_ADI_step(amiy, apiy, alphy, xs2y, rhs)
        
        # convergence criterium
        CONV = NSE_norm_L2(rhs, dx, dy)
        
        # computes u^{n+1}
        @views u[:] = u[:] + rhs[:]
        
        @printf("It=%d   time=%5.3f ||u-uold||=%10.5e \n", Niter, temps, CONV)
    end

    fex = NSE_fexact(Lx, Ly, xx, yy)

    @printf("It=%d   time=%5.3f ||u-uold|| = %10.5e \n", Niter, temps, CONV)
    @printf("Norm ||Uex-Unum|| = %10.5e \n",NSE_norm_L2(fex-u, dx, dy))

end

@time main()
@time main()