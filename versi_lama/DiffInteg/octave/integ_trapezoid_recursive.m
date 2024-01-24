function Ih = integ_trapezoid_recursive( func, a, b, I2h, k )
%
% Recursive trapezoidal rule
%
% func: handle of function being integrated
% a, b: limits of integration
% I2h : integral with 2^(k-1) panels
% Ih : integral with 2^k

if k == 1
  fa = feval( func, a );
  fb = feval( func, b );
  Ih = (fa + fb)*(b - a)*0.5;
else
  N = 2^(k-2);
  h = (b-a)/N;
  x = a + 0.5*h;
  s = 0.0;
  for i = 1:N
    s = s + feval(func,x);
    x = x + h;
  end
  Ih = 0.5*( I2h + h*s );
end




