function root = secant( f, x, tol, NiterMax )

  if ~exists("tol", "local")
    tol = 1e-9
  end

  if ~exists("NiterMax", "local")
    NiterMax = 30
  end

  x_old  = 0.0
  fx_old = 0.0

  printf("\nRoot searching via secant method\n")
  printf("Initial guess: [%18.10f]\n", x)
  printf("Tolerance: %18.10e\n", tol)
  printf("secant will iterate up to %d maximum iterations.\n", NiterMax)

  printf("secant: Begin iteration\n")
  for iter = 1:NiterMax

    fx = f(x)
    
    printf("Secant: %5d (%18.10f,%18.10e)\n", iter, x, fx )
    if abs(fx) < tol
      printf("Secant: Convergence achieved\n")
      root = x
      return
    end

    if abs( x - x_old ) < %eps
      printf("Secant: abs(x - x_old) is less than tol eps\n")
      root = x
      return
    end

    if iter == 1
      SMALL = 1e-10
      fx_old = f(x+SMALL)
      dfx = abs( fx - fx_old ) / SMALL
    else
      dfx = ( fx - fx_old ) / ( x - x_old )
    end

    if abs(dfx) < %eps
      printf("ERROR: very small derivative\n")
      root = x
      return
    end

    x_old = x
    fx_old = fx
    
    dx = -fx/dfx
    x = x + dx
  
  end

  printf("Too many iterations")

endfunction
