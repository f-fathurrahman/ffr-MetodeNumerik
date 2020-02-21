function diffusion_1d_explicit( α, xf, tf, u0x, bx0, bxf, Nx, Nt )
    
    Δx = xf/(Nx-1)
    x = collect(range(0.0, stop=xf, length=Nx))
    
    Δt = tf/(Nt-1)
    t = collect(range(0.0, stop=tf, length=Nt))

    u = zeros(Float64, Nx, Nt)
    
    # Aplikasi syarat awal
    for i in 1:Nx
        u[i,1] = u0x(x[i])
    end
    
    # Syarat batas
    for k in 1:Nt
        u[1,k] = bx0(t[k])
        u[Nx,k] = bxf(t[k])
    end
    
    # The Fourier mesh number
    F = α*Δt/Δx^2
    
    if F >= 0.5
        @printf("diffusion_1d_explicit:\n")
        @printf("WARNING: F is greater than 0.5%f\n", F)
        @printf("WARNING: The solution is not guaranteed to be stable !!\n")
    else
        @printf("diffusion_1d_explicit:\n")
        @printf("INFO: F = %f >= 0.5\n", F)
        @printf("INFO: The solution should be stable\n")
    end

    for n in 1:Nt-1
        for i in 2:Nx-1
            u[i,n+1] = F*( u[i+1,n] + u[i-1,n] ) + (1 - 2*F)*u[i,n]
        end
    end
    
    return u, x, t

end