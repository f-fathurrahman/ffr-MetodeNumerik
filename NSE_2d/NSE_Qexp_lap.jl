# global dx dy Lx Ly;
# global nxm nym ;
# global ip im jp jm ic jc;

include("common.jl")
include("NSE_fsource.jl")
include("NSE_calc_lap.jl")

function main()

    # Input parameters
    Lx = 1.0
    Ly = 2.0
    Nx = 21
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

    display(xc); println()
    display(ym); println()

    #  Initialization
    u = zeros(Float64, Nxm, Nym)
    du = zeros(Float64, Nxm, Nym)

    # Time step
    dt = 0.5/( 1/dx^2 + 1/dy^2 )
    println("dt = ", dt)

    # Time loop
    CONV = 1.0
    Nitermax = 10 #10000

    Niter = 0
    Temps = 0.0

    #tcpu = cputime;  % to estimate the computational CPU time: initialization

    while( (CONV > 1e-6) && (Niter <= Nitermax) )
    
        Niter = Niter + 1
        Temps = Temps + dt
    
        # compute the vector u^{n+1}-u^n
        du = dt*( NSE_fsource(Lx,Ly,xx,yy) + NSE_calc_lap(u, dx, dy, imm, ip, jp, jm, ic, jc) )
    
        # convergence criterium
        #CONV = NSE_norm_L2(du);
    
        ## solution u^{n+1}
        #u = u + du;
    
        #% check for convergence
        #if( mod(niter,10) == 0 ); 
        #    fprintf('It=%4d   time=%5.3f ||u-uold||=%10.5e \n', niter, temps, eps);
        #end;

        #convt(niter) = niter;
        #conve(niter) = eps;

    end
     
    #% Exact solution
    #fex = NSE_fexact(Lx,Ly,xx,yy);
         
    #fprintf('\n=====End of computation ======= CPU time =%d\n',cputime-tcpu)
    #fprintf('It=%d   time=%5.3f ||u-uold||=%10.5e \n', niter, temps, eps);
    #fprintf('Norm ||Uex-Unum|| =%10.5e \n', NSE_norm_L2(fex-u));

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

    println("Pass here ...")
end

main()