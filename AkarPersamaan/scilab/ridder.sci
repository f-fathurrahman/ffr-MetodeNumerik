function root = ridder( f, a, b, tol, NiterMax )

  if ~exists("tol", "local") then
    tol = 1e-9
  end

  if ~exists("NiterMax", "local") then
    NiterMax = 30
  end

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
    printf("ERROR: Root is not bracketed\n")
    return 0.0
  end

  xOld = 0.0

  for iter = 1:NiterMax

    c = 0.5*(a + b)
    fc = f(c)

    s = sqrt(fc^2 - fa*fb)

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

    printf("Ridder: %5d %18.10f %18.10e\n", iter, x, fx)

    // Test for Convergence
    if abs(fx) < tol
      printf("Ridder: Convergence achieved\n")
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
