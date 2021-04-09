using Printf

function main()
    Ndim = 2
    NelsX = 3
    NelsY = 4
    Nelements = NelsX*NelsY

    NnodesPerElement = 4
    NnodesX = NelsX + 1
    NnodesY = NelsY + 1
    NnodesTotal = NnodesX*NnodesY
    g_coord = zeros(Ndim,NnodesTotal)

    ip = 0
    dx = 0.5
    dy = 0.4
    for i in 1:NnodesX
        for j in 1:NnodesY
            ip = ip + 1
            g_coord[1,ip] = (i-1)*dx
            g_coord[2,ip] = (j-1)*dy
        end
    end
    println("Global coordinate:")
    display(g_coord); println()

    # Connectivity matrix
    g_num = zeros(Int64, NnodesPerElement,Nelements) 
    # grid of global node numbers
    gnumbers = reshape( collect(1:NnodesTotal), NnodesY, NnodesX )
    iel = 0
    # element counter
    for i in 1:NelsX
        for j in 1:NelsY
            iel = iel + 1
            # clockwise numbering
            g_num[1,iel] = gnumbers[j,i] # local node 1
            g_num[2,iel] = gnumbers[j+1,i] # local node 2
            g_num[3,iel] = gnumbers[j+1,i+1] # local node 3
            g_num[4,iel] = gnumbers[j,i+1] # local node 4
        end
    end
    println("Global node numbers")
    display(g_num); println()

    iel = 1
    num = g_num[:,iel]
    coord = g_coord[:,num]'
    println("num = ", num)
    println("coord = ", coord)

end

main()