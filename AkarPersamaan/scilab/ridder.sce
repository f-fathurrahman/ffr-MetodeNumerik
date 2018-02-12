function root = ridder( f, a, b, tol, NiterMax )

  if ~exists("tol", "local") then
    tol = 1e-9
  end

  if ~exists("NiterMax", "local") then
    NiterMax = 30
  end

  printf("\nRoot searching via Ridder''s method\n")
  printf("Initial search interval: [%18.10f,%18.10f]\n", a, b)
  printf("Tolerance: %18.10e\n", tol)
  printf("ridder will iterate up to %d maximum iterations.\n", NiterMax)

  fa = f(a)
  if abs(fa) <= tol
    root = a
    return
  end

  fb = f(b)
  if abs(fb) <= tol
    root = b
    return
  end

  if fa*fb > 0.0
    printf("ERROR in ridder: Root is not bracketed\n")
    return 0.0
  end

  xOld = 0.0

  printf("ridder: Begin iteration\n") 
  for iter = 1:NiterMax

    // midpoint
    c = 0.5*(a + b)
    fc = f(c)

    s = sqrt(fc^2 - fa*fb)

    // if s is zero we need to 
    if s==0
      printf("ERROR: s is zero\n")
      return 0.0
    end

    dx = (c-a)*fc/s
    if (fa-fb) < 0
      dx = -dx
    end

    x = c + dx
    fx = f(x)

    printf("ridder: %5d %18.10f %18.10e\n", iter, x, fx)

    // Test for Convergence
    if abs(fx) < tol
      printf("ridder: Convergence achieved\n")
      root = x
      return
    end

    xOld = x

    // Rebracket the root as tightly as possible
    if fc*fx > 0.0
      if fa*fx < 0.0
        b = x
        fb = fx
      else
        a = x
        fa = fx
      end
    else
      a = c
      b = x
      fa = fc
      fb = fx
    end

  end

  printf("Too many iterations")
  root = x

endfunction
