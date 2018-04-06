function [t,y] = ode_euler(f,tspan,y0,N)
// Euler's method to solve vector differential equation y'(t) = f(t,y(t))
// for tspan = [t0,tf] and with the initial value y0 and N time steps
  
  if (~exists("N", "local")) | (N <= 0)
    N = 100
  end
  
  if ~exists("tspan","local")
    y0 = 0
  end
  
  h = (tspan(2) - tspan(1))/N
  t = tspan(1) + [0:N]'*h
  
  // check y0, transpose if needed
  if size(y0,1) ~= 1
    if size(y0,2) > 1
      y0 = y0'
    end
  end

  Ndim = size(y0,2)
  y = zeros(N+1,Ndim)
  
  // initial value
  y(1,:) = y0

  for k = 1:N
    y(k+1,:) = y(k,:) + h*f(t(k),y(k,:))
  end

endfunction
