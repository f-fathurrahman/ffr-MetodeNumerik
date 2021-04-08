# Based on: Simpson - Practical Finite Element Modelling in Earth Science

using Printf
using LinearAlgebra
using SparseArrays
import PyPlot

const plt = PyPlot

function calc_analytic_solution!(κ, A, K, L, t, xgrid, Nterms, Texact)
    Npoints = size(xgrid,1)
    # Analytic solution
    sumv = zeros(Npoints)
    for i in 1:Npoints
        x = xgrid[i]
        sumv[i] = 0.0
        for ni in 0:Nterms
            et = exp( -κ*(2*ni+1)^2 * pi^2*t/4/L^2 )
            sumv[i] = sumv[i] + (-1)^ni/(2*ni+1)^3*cos((2*ni+1)/2/L*pi*x)*et
        end
        Texact[i] = A * L^2/(2*K) * (1.0 - x^2/L^2 - 32/pi^3 * sumv[i] )
    end
    return
end



function main()
    
    # -------------------
    # Physical parameters
    # -------------------
    seconds_per_yr = 60*60*24*365  # Number of seconds in one year (s)
    
    Lx = 1e4  # Length of spatial domain (m)

    Cp = 1e3  # rock heat capacity (J/Kg/K)

    rho = 2700  # rock density

    K = 3.3  # bulk thermal conductivity (W/m/K)

    κ = K/(Cp*rho)  # thermal diffusivity

    Tb = 0.0    # Temperature at boundaries (degree Celcius)

    A = 2.6e-6  # Heat production per volume per second (W/m^3)

    H = A/(rho*Cp)  # Heat source term (K/s)


    # --------------------
    # Numerical parameters
    # --------------------
    dt = 1000*seconds_per_yr
    Ntime = 5000
    Nelements  = 40  # Total number of elements
    NnodesPerElement = 2  # no. of nodes per element
    NnodesTotal = Nelements + 1  # Total number of nodes
    dx = Lx/Nelements

    # Spatial domain (uniform mesh)
    g_coord = 0.0:dx:Lx  # the size NnodesTotal

    # Boundary condition
    bcdof = [1 NnodesTotal] # boundary nodes
    bcval = [Tb Tb] # boundary values

    # Connectivity and equation numbering
    g_num = zeros(Int64,NnodesPerElement,Nelements);
    g_num[1,:] = collect(1:NnodesTotal-1)
    g_num[2,:] = collect(2:NnodesTotal)

    # Matrices and vectors
    ff = zeros(Float64, NnodesTotal, 1)  # system load vector
    b = zeros(Float64, NnodesTotal, 1)   # system rhs vector
    LHS = spzeros(Float64, NnodesTotal, NnodesTotal) # system lhs matrix
    RHS = spzeros(Float64, NnodesTotal, NnodesTotal) # system rhs matrix
    displ = zeros(Float64, NnodesTotal, 1)  # initial temperature

    # ---------------
    # Matrix assembly
    # ---------------

    for iel = 1:Nelements
    
        num = g_num[:,iel] # get equation number

        dx = abs( g_coord[num[2]] - g_coord[num[1]] ) # length of the element
    
        # Mass matrix
        MM = dx*[ 1/3 1/6;
                  1/6 1/3]

        # Stiffness
        KM = [ κ/dx -κ/dx;
              -κ/dx  κ/dx]

        # Load vector
        F = dx*H*[1/2; 1/2]

        LHS[num,num] = LHS[num,num] + MM/dt + KM

        RHS[num,num] = RHS[num,num] + MM/dt

        ff[num] = ff[num] + F
    end

    # For analytic solution
    Nterms = 100
    L = Lx/2
    NptsGrid = 100
    xgrid = range(-L, stop=L, length=NptsGrid)
    Texact = zeros(NptsGrid)

    Nii_plot = [1, 100, 500, 1000, 2500, 5000] # for plotting

    # ---------
    # Time loop
    # ---------
    t = 0.0
    k = 1
    
    LHS[1,2] = 0.0
    LHS[1,1] = 1.0
    LHS[NnodesTotal,2] = 0.0
    LHS[NnodesTotal,NnodesTotal] = 1.0
    
    displ = zeros(NnodesTotal,1)
    b = zeros(NnodesTotal,1)
    factorLHS = lu(LHS)

    for ni in 1:Ntime

        #@printf("ni = %d\n", ni)
        t = t + dt

        b[:] = RHS*displ + ff
        b[bcdof] .= bcval
    
        #displ = LHS\b
        ldiv!(displ, factorLHS, b)

        #if ( do_plot && (ni in Nii_plot) )
        #    calc_analytic_solution!(κ, A, K, L, t, xgrid, Nterms, Texact)
        #    plt.clf()
        #    plt.plot(g_coord, displ, marker="o", label="numerical")
        #    plt.plot(xgrid .+ L, Texact, label="exact") # shift exact solution by L
        #    title_plot = @sprintf("t = %.1f year", t/seconds_per_yr)
        #    file_plot = @sprintf("IMG_sol_%04d.png", ni)
        #    plt.ylim(0.0, 10.0)
        #    plt.title(title_plot)
        #    plt.savefig(file_plot, dpi=150)
        #end

    end

end

@time main()
@time main()
