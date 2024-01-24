function root = newton_raphson( f, df, a, tol, NiterMax )

  if ~exists("tol", "local")
    tol = 1e-9
  end

  if ~exists("NiterMax", "local")
    NiterMax = 30
  end

  printf("\nRoot searching via Newton-Raphson method\n")
  printf("Initial guess root: %18.10f\n", a)
  printf("Tolerance: %18.10e\n", tol)
  printf("newton_raphson will iterate up to %d maximum iterations.\n", NiterMax)

  x = a

  printf("\n")
  for iter = 1:NiterMax

    fx = f(x)

    printf("newton_raphson: %5d %18.10f %18.10e\n", iter, x, fx )
    if abs(fx) < tol
      printf("newton_raphson: Convergence achieved\n")
      root = x
      return
    end

    dfx = df(x)

    if abs(dfx) < %eps
      printf("ERROR in newton_raphson: very small derivative\n")
      root = x
      return
    end

    dx = -fx/dfx
    x = x + dx

  end

  printf("WARNING in newton_rapshon: Iterations does not converge\n")
endfunction
