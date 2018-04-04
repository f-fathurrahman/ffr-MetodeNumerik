function [t,y] = ode_ABM(f,tspan,y0,N,KC)
// Adams-Bashforth-Moulton method to solve vector d.e. y’(t) = f(t,y(t))
//  for tspan = [t0,tf] and with the initial value y0 and N time steps
//  using the modifier based on the error estimate depending on KC = 1/0

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
  
  for k = 1:4
    F(k,:) = f(t(k),y(k,:))
  end
  
  p = y(4,:)
  c = y(4,:)
  KC22 = KC*251/270
  KC12 = KC*19/270
  h24 = h/24
  h241 = h24*[1 -5 19 9]
  h249 = h24*[-9 37 -59 55]
  
  for k = 4:N
    p1 = y(k,:) + h249*F
    m1 = p1 + KC22*(c-p)
    c1 = y(k,:) + h241*[F(2:4,:); f(t(k + 1),m1)]
    y(k+1,:) = c1 - KC12*(c1 - p1)
    p = p1
    c = c1
    F = [F(2:4,:); f(t(k + 1),y(k + 1,:))]
  end

endfunction