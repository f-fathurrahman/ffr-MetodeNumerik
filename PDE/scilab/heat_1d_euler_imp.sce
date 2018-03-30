function [u,x,t] = heat_1d_euler_imp(a,xf,T,initialTemp,bx0,bxf,Nx,Nt)

// solve 1d heat equation using implicit Euler method
//
//   a u_xx = u_t
//
// for: 0 <= x <= xf
//      0 <= t <= T
//
// Initial condition: u(x,0) = it0(x)
//
// Boundary condition: u(0,t) = bx0(t), u(xf,t) = bxf(t)
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

  for i = 1:Nx+1
    u(i,1) = initialTemp(x(i))
  end

  for it = 1:Nt+1
    u([1 Nx+1],it) = [bx0(t(it)); bxf(t(it))]
  end

  r  = a*dt/dx/dx
  r2 = 1 + 2*r
  
  for i = 1:Nx-1
    A(i,i) = r2
    if i > 1
      A(i-1,i) = -r
      A(i,i-1) = -r
    end
  end

  for k=2:Nt+1
    b = [r*u(1,k); zeros(Nx-3,1);
    r*u(Nx+1,k)] + u(2:Nx,k-1);
    //u(2:Nx,k) = trid(A,b)
    u(2:Nx,k) = A\b
  end

endfunction
