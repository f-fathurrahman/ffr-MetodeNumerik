function diffusion_1d_CN(
    L::Float64, Nx::Int64, T::Float64, Nt::Int64,
    α::Float64, u0x, bx0, bxf, f  # f: source term
)

    Δx = L/(Nx-1)
    x = collect(range(0.0, stop=L, length=Nx))

    Δt = T/(Nt-1)
    t = collect(range(0.0, stop=T, length=Nt))

    u = zeros(Float64, Nx, Nt)

    for i in 1:Nx
        u[i,1] = u0x(x[i])
    end

    for k in 1:Nt
        u[1,k] = bx0(t[k])
        u[Nx,k] = bxf(t[k])
    end

    F = α*Δt/Δx^2

    # Build matrix A and vector b
    A = zeros(Float64, Nx, Nx)
    b = zeros(Float64, Nx)
    for i in 2:Nx-1
        A[i,i] = 1 + F
        A[i,i-1] = -0.5*F
        A[i,i+1] = -0.5*F
    end
    A[1,1] = 1.0
    A[Nx,Nx] = 1.0

    for n in 1:Nt-1
        for i in 2:Nx-1
            b[i] = u[i,n] + 0.5*F*( u[i-1,n] - 2*u[i,n] + u[i+1,n] ) + 0.5*( f(x[i],t[n]) + f(x[i],t[n+1]) )*Δt
        end
        b[1] = 0.0
        b[Nx] = 0.0
        # Solve the linear equations
        u[:,n+1] = A\b
    end
    return u, x, t

end