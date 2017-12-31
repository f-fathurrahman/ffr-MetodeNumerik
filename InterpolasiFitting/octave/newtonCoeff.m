function a = newtonCoeff(xData,yData)
% Returns coefficients of Newton's polynomial.
% USAGE: a = newtonCoeff(xData,yData)
% xData = x-coordinates of data points.
% yData = y-coordinates of data points.
n = length(xData);
a = yData;
for k = 2:n
  a(k:n) = (a(k:n) - a(k-1))./(xData(k:n) - xData(k-1));
end