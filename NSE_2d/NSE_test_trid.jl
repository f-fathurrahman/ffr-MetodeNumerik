using LinearAlgebra
using Random

include("NSE_trid_per_c2D.jl")

function main()

    #Random.seed!(1234)

    N = 10 # Dimension of the matrix of the system

    a = rand(Float64, N) # -1 subdiagonal
    c = rand(Float64, N) # +1 subdiagonal
    b = -(a + c) # 0 diagonal

    # build the matrix of the system
    # lower diagonal, diagonal, upper diagonal
    A = diagm( -1 => a[2:N], 0 => b, 1 => c[1:N-1] )
    A[1,N] = a[1]
    A[N,1] = c[N]
    
    # set the matrix diagonal dominant => invertible
    maxA = norm(A, Inf)
    println("maxA = ", maxA)

    display(A)
    println()

    b[:] = b[:] .+ maxA
    A = A + diagm( 0 => ones(N)*maxA )

    display(A)
    println()

    # manufactured solution
    sol_vec = collect(range(1.0, stop=N))
    
    # corrsponding RHS term
    fi = A*sol_vec
    
    display(fi)
    println()

    # Using \
    #sol_vec2 = A \ fi
    #display(sol_vec2)
    #println()

    # Solution using trid_per_2
    sol_trid = NSE_trid_per_c2D(a, b, c, fi)
    display(sol_trid)
    println()

    # Compare the solutions
    #for j=1:m
    #    fprintf('===== System  number  %i\n',j)
    #    fprintf('NSE_trid_per_c2D //   Matlab //   Exact\n')
    #    disp([soltrid(j,:)' solmat(:,j) sol(:,j)])
    #    fprintf('Press return to continue\n');pause
    #end

    println("Pass here")
end

main()
