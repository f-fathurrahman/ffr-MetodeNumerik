function [u,x,t] = wave_1d(a,xf,tf,it0,i1t0,bx0,bxf,Nx,Nt)
// solve a u_xx = u_tt for 0<=x<=xf, 0<=t<=T
//  Initial Condition: u(x,0) = it0(x), u_t(x,0) = i1t0(x)
//  Boundary Condition: u(0,t)= bx0(t), u(xf,t) = bxf(t)
//  M = # of subintervals along x axis
//  N = # of subintervals along t axis
  dx = xf/Nx
  x = [0:Nx]'*dx
  
  dt = tf/Nt
  t = [0:Nt]*dt
  
  r = a*(dt/dx)^2
  r1 = r/2
  r2 = 2*(1 - r)
  if r > 1
    printf("WARNING: propagation will not be stable\n\n")
    printf("r = %f > 1\n", r)
  end

  // initial conditions
  for i = 1:Nx+1
    u(i,1) = it0(x(i))
  end

  // boundary conditions
  for k = 1:Nt+1
    u(1,k)    = bx0(t(k))
    u(Nx+1,k) = bxf(t(k))
  end

  u(2:Nx,2) = r1*u(1:Nx-1,1) + (1-r)*u(2:Nx,1) + r1*u(3:Nx+1,1) + dt*i1t0(x(2:Nx))
  
  for k = 3:Nt+1
    u(2:Nx,k) = r*u(1:Nx-1,k-1) + r2*u(2:Nx,k-1) + r*u(3:Nx+1,k-1) - u(2:Nx,k-2)
  end

endfunction
