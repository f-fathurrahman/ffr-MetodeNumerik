using Printf

# Using column-first

function main()
    Ndim = 2
    Nx = 4
    Ny = 5
    Nels = Nx*Ny
    g_coord = zeros(Ndim,Nels)
    ip = 0
    dx = 0.5
    dy = 0.4
    for j in 1:Ny
        for i in 1:Nx
            ip = ip + 1
            g_coord[1,ip] = (i-1)*dx
            g_coord[2,ip] = (j-1)*dy
        end
    end
    display(g_coord); println()
    println("Pass here")
end

main()