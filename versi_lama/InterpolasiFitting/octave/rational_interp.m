function yInterp = rational_interp(xData,yData,x)
% Rational function interpolation;
% returns value of the interpolant at x.
% USAGE: yInterp = rational(xData,yData,x)
% xData = x-coordinates of data points.
% yData = y-coordinates of data points.
m = length(xData);
r = yData;
rOld = zeros(1,m);

for k = 1:m-1
  for i = 1:m-k
    if x == xData(i+k)
      yInterp = yData(i+k);
      return
    else
      c1 = r(i+1) - r(i);
      c2 = r(i+1) - rOld(i+1);
      c3 = (x - xData(i))/(x - xData(i+k));
      r(i) = r(i+1) + c1/(c3*(1.0 - c1/c2) - 1.0);
      rOld(i+1) = r(i+1);
    end
  end
end
yInterp = r(1);