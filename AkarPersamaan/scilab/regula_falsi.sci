function root = regula_falsi( f, x1, x2, tol, NiterMax )

  if ~exists("tol", "local") then
    tol = 1e-9
  end

  if ~exists("NiterMax", "local") then
    NiterMax = 100
  end

  f1 = f(x1)
  f2 = f(x2)

  if(f1*f2 > 0.0)
    printf("ERROR: Root is not bracketed")
    return
  end

  for iter = 1:NiterMax

    x3 = (x1*f2 - x2*f1)/(f2 - f1)
    f3 = f(x3)

    printf("Regula falsi: %5d (%18.10f,%18.10e)\n", iter, x3,f3 )

    if abs(f3) < tol
      printf("Regula falsi: Convergence achieved\n")
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

  root = x3

endfunction
