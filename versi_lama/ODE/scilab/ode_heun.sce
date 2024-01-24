function [t,y] = ode_heun(f,tspan,y0,N)
  // Heun method to solve vector differential equation y’(t) = f(t,y(t))
  // for tspan = [t0,tf] and with the initial value y0 and N time steps

  if (~exists("N", "local")) | (N <= 0)
    N = 100
  end
  
  if ~exists("tspan","local")
    y0 = 0
  end

  h = (tspan(2) - tspan(1))/N

  t = tspan(1) + [0:N]'*h

  y(1,:) = y0(:)' // initial value as row vector

  for k = 1:N
    fk = f(t(k),y(k,:))
    y(k+1,:) = y(k,:) + h*fk
    y(k+1,:) = y(k,:) + h/2*( fk + f(t(k+1),y(k+1,:)) )
  end

endfunction

