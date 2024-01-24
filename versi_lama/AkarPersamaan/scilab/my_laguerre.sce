function x = my_laguerre(a,tol)

  // start from random number
  x = rand()

  n = length(a) - 1
  MaxIter = 30

  for i = 1:MaxIter
    [p, dp, d_dp] = evalpoly(a, x)
    if abs(p) < tol
      return
    end
    g = dp/p
    h = g*g - d_dp/p
    f = sqrt( (n-1)*(n*h - g*g) )
    if abs(g + f) >= abs(g - f)
      dx = n/(g + f)
    else
      dx = n/(g - f)
    end
    x = x - dx
    if abs(dx) < tol
      return
    end
  end
  error("Too many iterations in my_laguerre")

endfunction

