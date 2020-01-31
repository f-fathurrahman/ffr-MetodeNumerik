function heat_1d_CN( α, xf, tf, u0x, bx0, bxf, Nx, Nt )

    dx = xf/Nx
    x = collect(range(0.0, stop=xf, length=Nx+1))

    dt = tf/Nt
    t = collect(range(0.0, stop=tf, length=Nt+1))

    u = zeros( Float64, Nx+1, Nt+1 )

    # Aplikasi syarat awal
    for i in 1:Nx+1
        u[i,1] = u0x( x[i] )
    end

    # Syarat batas
    for k in 1:Nt+1
        u[1,k] = bx0( t[k] )
        u[Nx+1,k] = bxf( t[k] )
    end

    r = α*dt/dx^2

    A = zeros( Float64, Nx-1, Nx-1 )
    for i in 1:Nx-1
        A[i,i] = 2*(1 + r)
        if i > 1
            A[i-1,i] = -r
            A[i,i-1] = -r
        end
    end

    B = zeros( Float64, Nx-1, Nx-1 )
    for i in 1:Nx-1
        B[i,i] = 2*(1 - r)
        if i > 1
            B[i-1,i] = r
            B[i,i-1] = r
        end
    end

    b = zeros(Float64, Nx-2)
    for k in 2:Nt+1
        b = B*u[2:Nx,k-1]
        u[2:Nx,k] = A\b
    end

    return u, x, t

end