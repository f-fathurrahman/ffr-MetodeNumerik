function sigma = stdDev(coeff,xData,yData)
% Returns the standard deviation between data points and the polynomial
% a(1)*x^(m-1) + a(2)*x^(m-2) + ... + a(m)
% USAGE: sigma = stdDev(coeff,xData,yData)
% coeff = coefficients of the polynomial.
% xData = x-coordinates of data points.
% yData = y-coordinates of data points.

m = length(coeff);
n = length(xData);

sigma = 0;

for i =1:n
  y = polyEval(coeff,xData(i));
  sigma = sigma + (yData(i) - y)^2;
end
sigma =sqrt(sigma/(n - m));

endfunction

function y = polyEval(coeff,x)
% Returns the value of the polynomial at x.
m = length(coeff);
y = coeff(1);
for j = 1:m-1
  y = y*x + coeff(j+1);
end
endfunction