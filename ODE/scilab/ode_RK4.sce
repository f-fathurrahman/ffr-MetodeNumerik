function [t,y] = ode_RK4(f,tspan,y0,N)

  if (~exists("N", "local")) | (N <= 0)
    N = 100
  end
  
  if ~exists("tspan","local")
    y0 = 0
  end

  // check y0, transpose if needed
  if size(y0,1) ~= 1
    if size(y0,2) > 1
      y0 = y0'
    end
  end

  Ndim = size(y0,2)
  y = zeros(N+1,Ndim)

  y(1,:) = y0
  h = (tspan(2) - tspan(1))/N
  t = tspan(1) + [0:N]'*h
  h2 = h/2
  
  for k=1:N
    f1 = h*f( t(k),y(k,:) )
    f2 = h*f( t(k) + h2, y(k,:) + f1/2)
    f3 = h*f( t(k) + h2, y(k,:) + f2/2)
    f4 = h*f( t(k) + h, y(k,:) + f3)
    y(k+1,:) = y(k,:) + (f1 + 2*(f2+f3) + f4)/6
  end

endfunction

