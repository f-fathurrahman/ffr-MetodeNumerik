function [u,x,t] = heat_1d_euler_exp(a,xf,T,initialT,bx0,bxf,Nx,Nt)

// solve 1d heat equation:
//
//   a u_xx = u_t
//
// for: 0<=x<=xf, 0<=t<=T
//
// Initial condition: u(x,0)=it0(x)
//
// Boundary condition: u(0,t)=bx0(t), u(xf,t)=bxf(t)
//
// Nx = no. of subintervals along x-axis
//      Total number of points is Nx + 1
//
// Nt = no. of subintervals along t-axis
//      Total number of time is Nt + 1

  dx = xf/Nx
  x  = [0:Nx]'*dx

  dt = T/Nt
  t  = [0:Nt]*dt

  // Set initial condition
  for i = 1:Nx+1
    u(i,1) = initialT(x(i))
  end
  
  for it = 1:Nt+1
    u([1 Nx+1],it) = [bx0(t(it)); bxf(t(it))]
  end
  
  r  = a*dt/dx/dx
  r1 = 1 - 2*r

  if r > 0.5
    printf("\nheat_1d_euler:\n")
    printf("WARNING: r is larger than 0.5: %f\n", r)
    printf("WARNING: solution is not stable\n\n")
  else
    printf("\nheat_1d_euler:\n")
    printf("r = %f >= 0.5\n", r)
    printf("The solution should be stable\n\n")
  end

  for it = 1:Nt
    for i = 2:Nx
      u(i,it+1) = r*( u(i+1,it) + u(i-1,it) ) + r1*u(i,it)
    end
  end
 
endfunction


