function [t,y]=ode_RK4(f,tspan,y0,N,varargin)
// Runge-Kutta method to solve vector differential equation y'(t)=f(t,y(t))
// for tspan=[t0,tf] and with the initial value y0 and N time steps

  if (~exists("N", "local")) | (N <= 0)
    N = 100
  end
  
  if ~exists("tspan","local")
    y0 = 0
  end

  y0 = y0(:)' // make it a row vector
  y = zeros(N+1,length(y0))
  y(1,:) = y0
  h = (tspan(2) - tspan(1))/N
  t = tspan(1) + [0:N]'*h
  h2 = h/2
  for k=1:N
    f1 = h*f(t(k),y(k,:))
    f1 = f1(:)
    f2 = h*f(t(k)+h2,y(k,:)+f1/2)
    f2 = f2(:)'
    f3 = h*f(t(k)+h2,y(k,:)+f2/2)
    f3 = f3(:)'
    f4 = h*f(t(k)+h,y(k,:)+f3)
    f4 = f4(:)'
    y(k+1,:) = y(k,:) + (f1 + 2*(f2+f3) + f4)/6
  end

endfunction

