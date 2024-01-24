% finding all roots within a given interval

func = @(x) (x - tan(x));

a = 0.0;
b = 20.0;
dx = 0.01;

Nroots = 0;

while true
  [x1, x2] = rootsearch(func, a, b, dx);
  if isnan(x1)
    break
  else
    a = x2; % next leftmost point for interval to be searched for roots
    x = bisection( func, x1, x2, true );
    if ~isnan(x)
      Nroots = Nroots + 1;
      % append new root
      root(Nroots) = x;
    end
  end
end
root