function u = poisson3d(u0, x, y, z, Nx, Ny, Nz, TOL, func)
// func is given as 3d-array
// Nx, Ny: no of nodes in x and y direction
  
  MaxIter = 10000

  hx = (x(Nx) - x(1))/(Nx-1)
  kx = 1.0/(hx*hx)

  hy = (y(Ny) - y(1))/(Nx-1)
  ky = 1.0/(hy*hy)

  hz = (z(Ny) - z(1))/(Nz-1)
  kz = 1.0/(hz*hz)
  
  kxyz = 2.0*(kx + ky + kz)

  u = u0
  for iter =1:MaxIter
    err = 0.0
    u_old = u
    // loop only for internal nodes
    for k = 2:Nz-1
    for j = 2:Ny-1
    for i = 2:Nx-1
        u(i,j,k) = ( kx*( u(i-1,j,k) + u(i+1,j,k) ) + ...
                     ky*( u(i,j-1,k) + u(i,j+1,k) ) + ...
                     kz*( u(i,j,k-1) + u(i,j,k+1) ) - func(i,j,k) ) / kxyz
    end
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
