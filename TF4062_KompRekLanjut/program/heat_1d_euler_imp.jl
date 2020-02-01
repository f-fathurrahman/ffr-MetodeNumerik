function heat_1d_euler_imp( α, xf, tf, u0x, bx0, bxf, Nx, Nt )

    dx = xf/(Nx-1)
    x = collect(range(0.0, stop=xf, length=Nx))

    dt = tf/(Nt-1)
    t = collect(range(0.0, stop=tf, length=Nt))

    u = zeros( Float64, Nx, Nt )

    # Aplikasi syarat awal
    for i in 1:Nx
        u[i,1] = u0x( x[i] )
    end

    # Syarat batas
    for k in 1:Nt
        u[1,k] = bx0( t[k] )
        u[Nx,k] = bxf( t[k] )
    end

    r = α*dt/dx^2

    # Bangun matriks A
    A = zeros( Float64, Nx-2, Nx-2 )
    for i in 1:Nx-2
        A[i,i] = 1 + 2*r
        if i > 1
            A[i-1,i] = -r
            A[i,i-1] = -r
        end
    end

    # Bangun vektor b
    b = zeros(Float64, Nx-2)
    for k in 2:Nt
        
        b[:] = u[2:Nx-1,k-1] # copy
        b[1] = b[1] + r*u[1,k]
        b[Nx-2] = b[Nx-2] + r*u[Nx,k]
        
        # Selesaikan sistem persamaan linear
        u[2:Nx-1,k] = A\b
    
    end
    return u, x, t

end