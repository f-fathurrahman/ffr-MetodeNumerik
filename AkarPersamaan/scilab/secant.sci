function root = secant( f, x, tol, NiterMax )

  if ~exists("tol", "local")
    tol = 1e-9
  end

  if ~exists("NiterMax", "local")
    NiterMax = 30
  end

  x_old  = 0.0
  fx_old = 0.0

  for iter = 1:NiterMax

    fx = f(x)

    printf("Secant: %5d (%18.10f,%18.10e)\n", iter, x, fx )
    if abs(fx) < tol
      printf("Secant: Convergence achieved\n")
      root = x
      return
    end

    dfx = ( fx - fx_old ) / ( x - x_old )

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
