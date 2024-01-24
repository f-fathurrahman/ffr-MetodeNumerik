//> Find root of equation $f(x) = 0$ using regula falsi.
//> \texttt{x1} and \texttt{x2} is initial search interval
//> \texttt{tol} is optional and is set to $10^{-9}$ by default.
function root = regula_falsi( f, x1, x2, tol, NiterMax )

  if ~exists("tol", "local") then
    tol = 1e-9
  end

  if ~exists("NiterMax", "local") then
    NiterMax = 100
  end

  printf("\nRoot searching via regula falsi method\n")
  printf("Initial search interval: [%18.10f,%18.10f]\n", x1, x2)
  printf("Tolerance: %18.10e\n", tol)
  printf("regula_falsi will iterate up to %d maximum iterations.\n", NiterMax)

  f1 = f(x1)
  f2 = f(x2)

  if abs(f1) <= tol
    root = x1
    return
  end

  if abs(f2) <= tol
    root = x2
    return
  end

  if(f1*f2 > 0.0)
    printf("ERROR: Root is not bracketed")
    return
  end
  
//>
//> Regula falsi iterations starts here.
//>

  printf("regula_falsi: Begin iteration\n")
  for iter = 1:NiterMax

    x3 = (x1*f2 - x2*f1)/(f2 - f1)
    f3 = f(x3)

    printf("regula_falsi: %5d %18.10f %18.10e\n", iter, x3,f3 )

    if abs(f3) < tol
      printf("regula_falsi: Convergence achieved in %d iterations\n", iter)
      root = x3
      return
    end

    if f2*f3 < 0.0
      x1 = x3
      f1 = f3
    else
      x2 = x3
      f2 = f3
    end

  end

endfunction
