using Printf

function main()
    Ndim = 2
    NnodesX = 4
    NnodesY = 5
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
    display(g_coord); println()
    println("Pass here")
end

main()