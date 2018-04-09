function [x,y,idx_div] = sch_numerov(E,V,xspan,y0,dy0,N)
//
// the differential equation is defined by scalar E and function V
//
// y'' = 2*(V(x) - E)
// y(0) = y0
// y'(0) = y'(0)
//
// div_idx is index of divergence point
  
  h = (xspan(2) - xspan(1))/N
  
  x = xspan(1) + [0:N]'*h // column vector
  y = zeros(N+1,1)  // column vector

  y(1) = y0
  dy = dy0

  // Euler-Cromer step to calculate y(2)
  d2y = 2.0*( V(x(1)) - E )*y(1)
  dy = dy + h*d2y
  y(2) = y(1) + h*dy

  h6 = h*h/6.0
  ux1 = 1 - h6*( V(x(1)) - E )
  ux  = 1 - h6*( V(x(2)) - E )

  // default return value for div_idx
  idx_div = N

  for i = 2:N
    up1 = 1.0 - h6*( V(x(i+1)) - E )
    y(i+1) = ( ( 12.0 - 10.0*ux )*y(i) - ux1*y(i-1) ) / up1
    ux1 = ux
    ux = up1
    //
    if abs(y(i+1)) > 1e10
      idx_div = i
      break
    end
  end

endfunction