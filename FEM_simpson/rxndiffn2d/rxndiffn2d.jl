using Printf
using LinearAlgebra
using SparseArrays

import PyPlot
const plt = PyPlot

import Random

function main()

    Random.seed!(1234)

    # physical parameters
    d = 20    # diffusivity of species B
    a = 0.05  # growth rate of species A
    b = 1.0   # growth rate of species B
    γ = 600    # kinetics
    NOISE_AMP = 0.01 # max. amplitude of random noise

    Lx = 5.0   # length of x domain
    Ly = 5.0   # length of y domain

    Ndim = 2
    NelsX = 50
    NelsY = 50
    Nelements = NelsX*NelsY
    #
    # Rectangular quadratic element
    NnodesPerElement = 9
    NnodesX = 2*NelsX + 1
    NnodesY = 2*NelsY + 1
    NnodesTotal = NnodesX*NnodesY

    dx = Lx/NelsX   # element length in x-direction
    dy = Ly/NelsY   # element length in y-direction

    NdofsPDE = 2  # number of degrees of freedom in pde
    NdofsPerElement = NdofsPDE*NnodesPerElement  # total degrees of freedom in one element

    # generate mesh, node and equation numbering
    
    # define mesh (numbering in y direction first)
    g_coord = zeros(Ndim,NnodesTotal)
    ip = 0
    for i in 1:NnodesX       # loop over nodes in x-direction
        for j in 1:NnodesY   # loop over nodes in y-direction
            ip = ip + 1
            g_coord[1,ip] = 0.5*(i-1)*dx
            g_coord[2,ip] = 0.5*(j-1)*dy
        end
    end
    #display(g_coord); println()

    # establish node numbering for each element
    # grid of global node numbers
    gnumbers = reshape( collect(1:NnodesTotal), NnodesY, NnodesX )
    # Connectivity matrix
    g_num = zeros(Int64, NnodesPerElement,Nelements)
    iel = 0
    for i in 1:2:NnodesX-1
        for j in 1:2:NnodesY-1
            iel = iel + 1
            g_num[1,iel] = gnumbers[j,i]
            g_num[2,iel] = gnumbers[j+1,i]
            g_num[3,iel] = gnumbers[j+2,i]
            g_num[4,iel] = gnumbers[j+2,i+1]
            g_num[5,iel] = gnumbers[j+2,i+2]
            g_num[6,iel] = gnumbers[j+1,i+2]
            g_num[7,iel] = gnumbers[j,i+2]
            g_num[8,iel] = gnumbers[j,i+1]
            g_num[9,iel] = gnumbers[j+1,i+1]
        end
    end
    #display(g_num); println()

    # establish equation number for each node
    idof = 0
    # system degrees of freedom (dof) counter
    nf = zeros(Int64, NdofsPDE, NnodesTotal) # node degree of freedom array
    for n in 1:NnodesTotal
        for i in 1:NdofsPDE
            idof = idof + 1
            nf[i,n] = idof  # store eqn number on each node
        end
    end
    NdofsTotal = idof
    
    # equation number for each element
    g = zeros(Int64, NdofsPerElement)   # equation numbers for 1 element
    g_g = zeros(Int64, NdofsPerElement, Nelements) # equation numbers for all elements
    # loop over elements
    for iel in 1:Nelements
        num = g_num[:,iel]  # node numbers for this element
        # initialise local eqn number
        # loop 2 times over nodes of an element, once for each dof
        inc = 0
        for i in 1:NnodesPerElement
            inc = inc + 1
            g[inc] = nf[1,num[i]]
        end
        for i in 1:NnodesPerElement
            inc = inc+1
            g[inc] = nf[2,num[i]]
        end
        g_g[:,iel] = g  # store the equation numbers
    end
    #display(g_g); println()

    #
    # FIXME: no boundary conditions?
    #


    NintegPoints  = 9 # number of integration points in an element
    integPoints = zeros(NintegPoints,Ndim)
    integPoints[1:3:7,1] .= -sqrt(0.6)
    integPoints[2:3:8,1] .= 0.0
    integPoints[3:3:9,1] .= sqrt(0.6)
    integPoints[1:3,2] .= sqrt(0.6)
    integPoints[4:6,2] .= 0.0
    integPoints[7:9,2] .= -sqrt(0.6)

    # Gauss weights
    w = [ 5/9, 8/9, 5/9]
    v = hcat( (5/9) .* w, (8/9) .* w, (5/9)*w )
    wIntegPoints = v[:]


    fun_s = zeros(NintegPoints, NintegPoints)
    der_s = zeros(Ndim, NintegPoints, NintegPoints)
    der = zeros(Ndim, NintegPoints)
    #
    # evaluate shape functions and their derivatives
    # at integration points and save the results
    #
    for k in 1:NintegPoints
        #
        ξ = integPoints[k,1]
        η = integPoints[k,2]
        ηm = η - 1
        ηp = η + 1
        ξm = ξ - 1
        ξp = ξ + 1
        x2p1 = 2*ξ + 1
        x2m1 = 2*ξ - 1
        e2p1 = 2*η + 1
        e2m1 = 2*η - 1
        # shape functions
        fun = [ 0.25*ξ*ξm*η*ηm, -0.5*ξ*ξm*ηp*ηm,
                0.25*ξ*ξm*η*ηp, -0.5*ξp*ξm*η*ηp,
                0.25*ξ*ξp*η*ηp, -0.5*ξ*ξp*ηp*ηm,
                0.25*ξ*ξp*η*ηm, -0.5*ξp*ξm*η*ηm,
                ξp*ξm*ηp*ηm ]
        # derivates of shape functions
        der[1,1] =  0.25*x2m1*η*ηm
        der[1,2] = -0.5*x2m1*ηp*ηm
        der[1,3] =  0.25*x2m1*η*ηp
        der[1,4] = -ξ*η*ηp
        der[1,5] =  0.25*x2p1*η*ηp
        der[1,6] = -0.5*x2p1*ηp*ηm
        der[1,7] =  0.25*x2p1*η*ηm
        der[1,8] = -ξ*η*ηm
        der[1,9] =  2*ξ*ηp*ηm
        der[2,1] = 0.25*ξ*ξm*e2m1
        der[2,2] = -ξ*ξm*η;
        der[2,3] = 0.25*ξ*ξm*e2p1 
        der[2,4] = -0.5*ξp*ξm*e2p1
        der[2,5] = 0.25*ξ*ξp*e2p1 
        der[2,6] = -ξ*ξp*η
        der[2,7] = 0.25*ξ*ξp*e2m1 
        der[2,8] = -0.5*ξp*ξm*e2m1
        der[2,9] = 2*ξp*ξm*η
        # save
        fun_s[k,:] = fun
        der_s[:,:,k] = der
    end

    LHS = spzeros(Float64, NdofsTotal, NdofsTotal)
    bv = zeros(Float64, NdofsTotal)
    displ = zeros(Float64, NdofsTotal)
    displ0 = zeros(Float64, NdofsTotal)
    displ_tmp = zeros(Float64, NdofsTotal)

    # initial conditions
    displ[nf[1,:]] = (a + b)*ones(NnodesTotal) + NOISE_AMP*randn(NnodesTotal)
    displ[nf[2,:]] = b/(a+b)^2 * ones(NnodesTotal) + NOISE_AMP*randn(NnodesTotal)

    ## dirichlet boundary conditions
    ## (zero-flux when the vectors are left empty)
    #bcdof = [ ] ; % fixed nodes
    #bcval = [ ] ; % fixed values
    # XXXX NOT USED?

    t = 0.0
    dt   = 5e-4  # time step
    Ntime = 1
    
    KMa = zeros(NnodesPerElement,NnodesPerElement)
    KMb = zeros(NnodesPerElement,NnodesPerElement)
    MM  = zeros(NnodesPerElement,NnodesPerElement)
    FA  = zeros(NnodesPerElement)
    FB  = zeros(NnodesPerElement)

    deriv = zeros(2,9)
    jac = zeros(2,2)
    invjac = zeros(2,2)

    for it in 1:200
        displ0[:] = displ[:] # save old solution
        t = t + dt  # update time
        @printf("t = %18.10f s\n", t)
        err = 1.0 # initialise error to arbitary large number
        #
        iterInner = 0
        while err > 1e-3
            iterInner = iterInner + 1
            #    
            bv[:] .= 0.0  # system rhs vector
            LHS.nzval[:] .= 0.0
            #
            # element integration and assembly
            #
            for iel in 1:Nelements
                num = g_num[:,iel] # node numbers
                ga = g_g[1:NnodesPerElement,iel] # equation numbers for A
                gb = g_g[NnodesPerElement+1:end,iel]          # equation numbers for B
                coord = g_coord[:,num]'  # node coordinates
                #
                fill!(KMa,0.0)
                fill!(KMb,0.0)
                fill!(MM,0.0)
                fill!(FA,0.0)
                fill!(FB,0.0)
                #
                for k in 1:NintegPoints 
                    fun = fun_s[k,:]  # shape functions
                    der = der_s[:,:,k] # derivs. of N in local coords
                    jac[:] = der*coord   # Jacobian matrix
                    detjac = det(jac)  # determinant of Jac
                    invjac[:] = inv(jac)  # inverse of Jac
                    #
                    deriv[:] = invjac*der # derivs. of N in physical coords
                    #
                    Ai = fun'*displ[ga]  # interpolate A to integration pt.
                    Bi = fun'*displ[gb]  # interpolate B to integration pt.
                    #
                    dwt = detjac*wIntegPoints[k] # multiplier
                    MM = MM + fun*fun'*dwt         # mass matrix
                    KMa = KMa + deriv'*deriv*dwt    # A diffn matrix
                    KMb = KMb + d*deriv'*deriv*dwt # B diffn matrix
                    
                    FA = FA + γ*(a + Ai^2*Bi)*fun*dwt # A load vector
                    FB = FB + γ*(b - Ai^2*Bi)*fun*dwt # B load vector
                end
                # assemble global lhs matrix and rhs vector
                for j in 1:9, i in 1:9
                    ia = ga[i]
                    ja = ga[j]
                    LHS[ia,ja] = LHS[ia,ja] + MM[i,j]/dt + KMa[i,j] + γ*MM[i,j]  # A contrib.
                end
                for j in 1:9, i in 1:9
                    ib = gb[i]
                    jb = gb[j]
                    LHS[ib,jb] = LHS[ib,jb] + MM[i,j]/dt + KMb[i,j]  # B contrib.
                end
                #LHS[ga,ga] = LHS[ga,ga] + MM/dt + KMa + γ*MM  # A contrib.
                #LHS[gb,gb] = LHS[gb,gb] + MM/dt + KMb         # B contribution
                bv[ga] = bv[ga] + MM*displ0[ga]/dt + FA  # A contribution
                bv[gb] = bv[gb] + MM*displ0[gb]/dt + FB  # B contribution
            end
            #println("Finish assembly")
            #%-------------------------------------------------
            #% implement boundary conditions and solve system
            #%-------------------------------------------------
            #% apply boundary conditions
            #lhs(bcdof,:) = 0;
            #tmp = spdiags(lhs,0) ;
            #tmp(bcdof)=1;
            #lhs=spdiags(tmp,0,lhs) ;
            #bv(bcdof) = bcval ;
            displ_tmp[:] = displ[:] # save solution vector
            displ = LHS\bv # solve sytem
            #% check for convergence
            err = maximum( abs.(displ - displ_tmp) )/maximum(abs.(displ))
            @printf("Inner iteration: %3d err = %18.5e\n", iterInner, err)
        end # end of nonlinear iteration loop
    end

    xg = reshape(g_coord[1,:],NnodesY,NnodesX)
    yg = reshape(g_coord[2,:],NnodesY,NnodesX)
    Ag = reshape(displ[nf[1,:]],NnodesY,NnodesX)
    Bg = reshape(displ[nf[2,:]],NnodesY,NnodesX)

    plt.clf()
    plt.pcolor(xg, yg, Ag, cmap="jet", shading="auto")
    plt.axis("equal")
    plt.title("A")
    plt.savefig("IMG_A.png", dpi=150)

    plt.clf()
    plt.pcolor(xg, yg, Bg, cmap="jet", shading="auto")
    plt.axis("equal")
    plt.title("B")
    plt.savefig("IMG_B.png", dpi=150)

    println("Pass here")
end

@time main()
