function root = newton_raphson( f, df, a, tol, NiterMax )

  if ~exists("tol", "local")
    tol = 1e-9
  end

  if ~exists("NiterMax", "local")
    NiterMax = 30
  end

  x = a

  for iter = 1:NiterMax

    fx = f(x)

    printf("Newton-Raphson: %5d (%18.10f,%18.10e)\n", iter, x, fx )
    if abs(fx) < tol
      printf("Newton-Raphson: Convergence achieved\n")
      root = x
      return
    end

    dfx = df(x)

    if abs(dfx) < %eps
      printf("ERROR: very small derivative\n")
      root = x
      return
    end

    dx = -fx/dfx
    x = x + dx

  end

  printf("Too many iterations")
endfunction
