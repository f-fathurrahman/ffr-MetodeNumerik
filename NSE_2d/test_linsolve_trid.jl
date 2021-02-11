using LinearAlgebra

include("linsolve_trid.jl")

function main()
    N = 5
    a =  1.0*ones(N)
    b = -2.0*ones(N)
    c =  1.0*ones(N)

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

    f = 10*ones(N)

    println("fullmat:")
    display(fullmat); println()

    println("inverse fullmat:")
    display(inv(fullmat)); println()

    x = zeros(N)
    linsolve_trid!(a, b, c, f, x)
    display(x); println()

    println("fullmat*x = ")
    println(fullmat*x)

    println("f = ")
    println(f)
end

main()