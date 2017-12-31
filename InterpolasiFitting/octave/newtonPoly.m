function p = newtonPoly(a,xData,x)
% Returns value of Newton's polynomial at x.
% USAGE: p = newtonPoly(a,xData,x)
% a = coefficient array of the polynomial;
%     must be computed first by newtonCoeff.
% xData = x-coordinates of data points.
n = length(xData);
p = a(n);
for k = 1:n-1;
  p = a(n-k) + (x - xData(n-k))*p;
end
