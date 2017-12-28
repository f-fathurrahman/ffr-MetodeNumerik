function root = bisection( f, x1, x2, tol )

  if ~exists("tol", "local") then
    tol = 1e-9
  end

  f1 = f(x1)
  f2 = f(x2)

  if(f1*f2 > 0.0)
    printf("ERROR: Root is not bracketed")
    return
  end

  Niter = int32( ceil( log(abs(x1-x2)/tol)/log(2.0) ) )
  printf("Bisection will iterate to %d iterations.\n", Niter)

  for iter = 1:Niter

    x3 = 0.5*(x1 + x2)
    f3 = f(x3)

    if( (abs(f3) > abs(f1)) & (abs(f3) > abs(f2) ) )
      printf("ERROR: x1, x2, x3: %f, %f, %f\n", x1, x3, x2)
      return
    end

    printf("Bisection: %5d (%18.10f,%18.10e)\n", iter, x3,f3 )

    if abs(f3) < tol
      printf("Bisection: Convergence achieved\n")
      break
    end

    if f2*f3 < 0.0
      x1 = x3
      f1 = f3
    else
      x2 = x3
      f2 = f3
    end

  end

  root = 0.5*(x1 + x2)
endfunction
