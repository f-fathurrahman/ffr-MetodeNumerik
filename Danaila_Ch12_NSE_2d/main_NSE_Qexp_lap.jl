using Printf
using LinearAlgebra: norm

include("common.jl")
include("NSE_calc_lap.jl")
include("NSE_norm_L2.jl")

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

    #display(xc); println()
    #display(ym); println()

    #  Initialization
    u = zeros(Float64, Nxm, Nym)
    du = zeros(Float64, Nxm, Nym)

    # Time step
    dt = 0.5/( 1/dx^2 + 1/dy^2 )
    println("dt = ", dt)

    # Time loop
    CONV = 1.0
    Nitermax = 100_000

    Niter = 0
    Temps = 0.0

    #tcpu = cputime;  % to estimate the computational CPU time: initialization

    convt = Int64[]
    conve = Float64[]

    while( (CONV > 1e-8) && (Niter <= Nitermax) )
    
        Niter = Niter + 1
        Temps = Temps + dt
    
        # compute the vector u^{n+1}-u^n
        du = dt*( NSE_fsource(Lx,Ly,xx,yy) + NSE_calc_lap(u, dx, dy, imm, ip, jp, jm, ic, jc) )
    
        # convergence criterium
        CONV = NSE_norm_L2(du, dx, dy)
    
        # solution u^{n+1}
        @views u[:,:] = u[:,:] + du[:,:]
    
        if( mod(Niter,10) == 0 ); 
            @printf("It=%4d  time=%5.3f norm(du) = %10.5e\n", Niter, Temps, CONV)
        end

        # Need this?
        append!(convt, Niter)

    end
     
    #% Exact solution
    fex = NSE_fexact(Lx, Ly, xx, yy)
    
    @printf("Norm ||Uex-Unum|| =%10.5e \n", NSE_norm_L2(fex - u, dx, dy))
    println( norm(fex-u)*sqrt(dx*dy) )

    #% Plot the iso-contours
    #% numerical sol./exact sol.
    #figure('Position',0.75*get(0,'Screensize'));
    #NSE_visu_isos(xx,yy,u,fex,11);
          
    #figure('Position',0.75*get(0,'Screensize'));
    #semilogy(convt,conve,'r-','LineWidth',2);
    #set(gca,'FontSize',24);
    #xlabel('Niter');
    #ylabel('\epsilon');
    #title('Convergence of the explicit solution','FontSize',24)

end

main()