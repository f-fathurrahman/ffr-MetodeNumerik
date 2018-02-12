// root of equation: x = f(x)
function root = fixed_point( f, x, tol, NiterMax )
  
  if ~exists("tol", "local") then
    tol = 1e-9
  end

  if ~exists("NiterMax", "local") then
    NiterMax = 100
  end

  for iter = 1:NiterMax
    fx = f(x)
    diffx = abs(x - fx)
    printf("fixed_point: %5d %18.10f %18.10f %18.10e\n", iter, x, fx, diffx)
    if diffx < tol
      printf("fixed_point: Convergence achieved in %d iterations\n", iter)      
      root = x
      return
    end
    // no convergence, prepare for next iteration
    x = fx
  end


endfunction