using Printf
using LinearAlgebra
using SparseArrays

function main()

    # physical parameters
    d = 20    # diffusivity of species B
    a = 0.05  # growth rate of species A
    b = 1.0   # growth rate of species B
    γ = 600    # kinetics
    NOISE_AMP = 0.01 # max. amplitude of random noise

    Lx = 5.0   # length of x domain
    Ly = 5.0   # length of y domain

    Ndim = 2
    NelsX = 2
    NelsY = 2
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


    println("Pass here")
end

@time main()