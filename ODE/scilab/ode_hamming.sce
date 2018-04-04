function [t,y] = ode_hamming(f,tspan,y0,N,KC)
// Hamming method to solve vector d.e. y’(t) = f(t,y(t))
// for tspan = [t0,tf] and with the initial value y0 and N time steps
// using the modifier based on the error estimate depending on KC = 1/0
  
  if ~exists("KC", "local")
    KC = 1
  end

  if ~exists("N","local")
    N = 100
  end
  
  if ~exists("tspan","local")
    y0 = 0
  end

  y0 = y0(:)'
  h = (tspan(2) - tspan(1))/N
  tspan0 = tspan(1) + [0 3]*h

  [t,y] = ode_RK4(f,tspan0,y0,3)
  t = [t(1:3)' t(4):h:tspan(2)]'
  for k = 2:4
    F(k-1,:) = f(t(k),y(k,:))
  end
  p = y(4,:)
  c = y(4,:)
  h34 = h/3*4
  KC11 = KC*112/121
  KC91 = KC*9/121
  h312 = 3*h*[-1 2 1]
  for k = 4:N
    p1 = y(k - 3,:) + h34*(2*(F(1,:) + F(3,:)) - F(2,:))
    m1 = p1 + KC11*(c-p)
    c1 = (-y(k-2,:) + 9*y(k,:) + h312*[F(2:3,:); f(t(k + 1),m1)])/8
    y(k+1,:) = c1 - KC91*(c1 - p1)
    p = p1
    c = c1
    F = [F(2:3,:); f(t(k+1),y(k+1,:))]
  end

endfunction