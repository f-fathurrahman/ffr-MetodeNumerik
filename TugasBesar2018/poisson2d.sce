function u = poisson2d(u0, x, y, Nx, Ny, TOL, func)
// Nx, Ny: no of nodes in x and y direction
  
  MaxIter = 10000

  hx = (x(Nx) - x(1))/(Nx-1)
  kx = 1.0/(hx*hx)

  hy = (y(Ny) - y(1))/(Nx-1)
  ky = 1.0/(hy*hy)

  kxy = 2.0*(kx + ky)

  // calculate array for RHS
  f = zeros(Nx,Ny)
  for j = 1:Ny
    for i = 1:Nx
      f(i,j) = func( x(i), y(j) )
    end
  end

  u = u0
  for iter =1:MaxIter
    err = 0.0
    u_old = u
    // loop only for internal nodes
    for j = 2:Ny-1
      for i = 2:Nx-1
        u(i,j) = ( kx*( u(i-1,j) + u(i+1,j) ) + ...
                   ky*( u(i,j-1) + u(i,j+1) ) - f(i,j) ) / kxy
      end
    end
    // calculate error
    err = sum(abs(u - u_old))
    printf("iter = %8d, err = %18.10e\n", iter, err)
    if err < TOL
      printf("Convergence is achieved\n")
      break
    end
  end

endfunction
