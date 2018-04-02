function [u,x,t] = heat_1d_CN(a,xf,T,initialTemp,bx0,bxf,Nx,Nt)

// solve 1d heat equation using explicit Crank-Nicholson
//
//   a u_xx = u_t
//
// for: 0<=x<=xf, 0<=t<=T
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
  r1 = 2*(1-r)
  r2 = 2*(1+r)

  // Build matrix A
  A = zeros(Nx-1,Nx-1)
  for i = 1:Nx-1
    A(i,i) = 2*(1 + r)
    if i > 1
      A(i-1,i) = -r
      A(i,i-1) = -r
    end
  end

  // Build matrix B
  B = zeros(Nx-1,Nx-1)
  for i = 1:Nx-1
    B(i,i) = 2*(1 - r)
    if i > 1
      B(i-1,i) = r
      B(i,i-1) = r
    end
  end

  for it = 2:Nt+1
    b = B*u(2:Nx,it-1)
    u(2:Nx,it) = A\b
  end

endfunction

