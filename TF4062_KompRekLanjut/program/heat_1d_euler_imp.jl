function heat_1d_euler_imp( α, xf, tf, u0x, bx0, bxf, Nx, Nt )

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

    # Bangun matriks A
    A = zeros( Float64, Nx-1, Nx-1 )
    for i in 1:Nx-1
        A[i,i] = 1 + 2*r
        if i > 1
            A[i-1,i] = -r
            A[i,i-1] = -r
        end
    end

    # Bangun vektor b
    b = zeros(Float64, Nx-1)
    for k in 2:Nt+1
        
        b[:] = u[2:Nx,k-1] # copy
        
        b[1] = b[1] + r*u[1,k]
        b[Nx-1] = b[Nx-1] + r*u[Nx+1,k]
        
        # Selesaikan sistem persamaan linear
        u[2:Nx,k] = A\b
    
    end
    return u, x, t

end