using LinearAlgebra

include("linsolve_trid.jl")

function main()
    N = 5
    a =  1.0*rand(Float64, N)
    b = -3.0*rand(Float64, N)
    c =  0.5*rand(Float64, N)

    # Build full matrix
    A = diagm(-1 => a[2:N], 0 => b, 1 => c[1:N-1])
    #
    println("A (full matrix):")
    display(A); println()

    # RHS vector
    f = 10*ones(Float64, N)

    println("inverse A:")
    display(inv(A)); println()

    x = zeros(Float64, N)
    linsolve_trid!(a, b, c, f, x)

    mae = sum(abs.(A*x - f))/length(f)
    println("Check solution, MAE (should be close to zero): ", mae)
end

main()