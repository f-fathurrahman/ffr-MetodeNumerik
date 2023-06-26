using LinearAlgebra

include("linsolve_trid.jl")
include("linsolve_trid_per.jl")

function main()
    N = 5
    a =  1.0*ones(N)
    b = -3.0*ones(N)
    c =  0.5*ones(N)

    fullmat = zeros(N,N)
    fullmat[1,1] = b[1]
    fullmat[1,2] = c[1]
    for i in 2:N-1
        fullmat[i,i-1] = a[i]
        fullmat[i,i] = b[i]
        fullmat[i,i+1] = c[i]
    end
    fullmat[N,N-1] = a[N]
    fullmat[N,N] = b[N]
    # periodic
    fullmat[N,1] = c[1]
    fullmat[1,N] = a[N]
    f = 10*ones(N)

    println("fullmat:")
    display(fullmat); println()

    println("inverse fullmat:")
    display(inv(fullmat)); println()

    x = zeros(N)
    linsolve_trid_per!(a, b, c, f, x)
    display(x); println()

    println("fullmat*x = ")
    println(fullmat*x)

    println("f = ")
    println(f)
end

main()