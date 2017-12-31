function yInterp = neville(xData,yData,x)
% Neville's polynomial interpolation;
% returns the value of the interpolant at x.
% USAGE: yInterp = neville(xData,yData,x)
% xData = x-coordinates of data points.
% yData = y-coordinates of data points.
n = length(xData);
y = yData;
for k = 1:n-1
  y(1:n-k) = ((x - xData(k+1:n)).*y(1:n-k) + (xData(1:n-k) - x).*y(2:n-k+1))./(xData(1:n-k) - xData(k+1:n));
end
yInterp = y(1);