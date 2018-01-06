//> Find root of equation $f(x) = 0$ using bisection.
//> \texttt{x1} and \texttt{x2} is initial search interval
//> \texttt{tol} is optional and is set to $10^{-9}$ by default.
function root = bisection( f, x1, x2, tol )

  if ~exists("tol", "local") then
    tol = 1e-9
  end

  Niter = int32( ceil( log(abs(x1-x2)/tol)/log(2.0) ) )
  
  printf("\nRoot searching via bisection algorithm\n")
  printf("Initial search interval: [%18.10f,%18.10f]\n", x1, x2)
  printf("Tolerance: %18.10e\n", tol)
  printf("bisection will iterate to %d iterations.\n", Niter)

//> Check whether the root is bracked or not.
//> The sign of \texttt{f1} and \texttt{f2} must be opposite.
  f1 = f(x1)
  f2 = f(x2)

  if(f1*f2 > 0.0)
    printf("ERROR in bisection: Root is not bracketed\n")
    return
  end

  printf("\n")

//>
//> Bisection iterations starts here.
//>
  for iter = 1:Niter

    // midpoint
    x3 = 0.5*(x1 + x2)
    // function value at the midpoint
    f3 = f(x3)

    // absolute value of f3 should not exceed f1 and f2 as it should
    if( (abs(f3) > abs(f1)) & (abs(f3) > abs(f2) ) )
      printf("ERROR: x1, x2, x3: %f, %f, %f\n", x1, x3, x2)
      return
    end

    printf("bisection: %5d %18.10f %18.10e\n", iter, x3,f3 )

    // check for convergence
    if abs(f3) < tol
      printf("bisection: Convergence achieved\n")
      break
    end

//>
//> Determine new search interval.
//>
    if f2*f3 < 0.0
      // root lies between f2 and f3
      // replace x1 with x2
      x1 = x3
      f1 = f3
    else
      // root lies between f1 and f3
      // replace x2 with x1
      x2 = x3
      f2 = f3
    end

  end

  root = 0.5*(x1 + x2)
endfunction
