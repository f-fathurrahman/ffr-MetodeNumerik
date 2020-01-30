function NSE_trid_per_c2D(a_, b_, c_, fi_)
    
    a = copy(a_)
    b = copy(b_)
    c = copy(c_)
    fi = copy(fi_)

    N = size(a, 1)
    
    # Modification of the matrix (1st and last coeff. on the diag.)
    b[1] = b[1] - a[1]
    b[N] = b[N] - c[N]

    # Initialization of the RHS for the calculation of  X2
	xs2 = zeros(Float64, N)
    xs2[1] = a[1]
    xs2[N] = c[N]

    # Resolution of the tridiagonal system
    # aa(j,i)*X(j,i-1) + ab(j,i)*X(j,i) + ac(j,i)*X(j,i+1) = xs2(j,i)
    #      i=1,n
    # -> solution stored in xs2

    #
    #--> Thomas algorithm 
    #-->                  alph(i) = 1/beta(i)
    #                     aa(i) = aa(i)/beta(i-1)
    #                     ac(i) = ac(i)/beta(i)

    alph = zeros(Float64, N)
    alph[1] = 1.0/b[1]
    for i in 2:N
        alph[i] = 1.0 / ( b[i] - a[i] * c[i-1] * alph[i-1] )
        a[i]    = a[i] * alph[i-1]
        c[i-1]  = c[i-1] * alph[i-1]
	end

    # -->resolution of the system [M*] X2 = v1_x
    #  ab(i)=beta(i)*gamma(i)
    b[1] = xs2[1]
    for i in 2:N
	   b[i] = xs2[i] - a[i] * b[i-1]
    end

    # xs2(i) = X2(i) 
    xs2[N] = b[N] * alph[N]
    
    for i in N-1:-1:1
        xs2[i] = b[i] * alph[i] - c[i] * xs2[i+1]
    end

    #display(fi)
    #println()

    # Resolution of the tridiagonal system
    # a(j,i)*X(j,i-1) + b(j,i)*X(j,i) + ac(j,i)*X(i+1)=fi(j,i), i=1,n
    # solution stored in fi
    # Forward substitution
    for i in 2:N
	   fi[i] = fi[i] - a[i] * fi[i-1]
    end
    #display(fi)
    #println()

    # Backward substitution
    fi[N] = fi[N] * alph[N]
    for i in N-1:-1:1
        fi[i] = fi[i] * alph[i] - c[i] * fi[i+1]
    end
    #display(fi)
    #println()

    # Final solution   	-->  v2_x = X*
    # dudx1 = X = [X1] - [X*] [X2]
    xs1 = ( fi[1] + fi[N] ) / ( 1.0 + xs2[1] + xs2[N] )

    for  i in 1:N
      fi[i] = fi[i] - xs1 * xs2[i]
    end
    #display(fi)
    #println()

    return fi
end
