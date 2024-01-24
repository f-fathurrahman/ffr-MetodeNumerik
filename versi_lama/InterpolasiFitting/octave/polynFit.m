function coeff = polynFit(xData,yData,m)
% Returns the coefficients of the polynomial
% a(1)*x^(m-1) + a(2)*x^(m-2) + ... + a(m)
% that fits the data points in the least squares sense.
% USAGE: coeff = polynFit(xData,yData,m)
% xData = x-coordinates of data points.
% yData = y-coordinates of data points.
A = zeros(m);
b = zeros(m,1);
s = zeros(2*m-1,1);

for i = 1:length(xData)
  temp = yData(i);
  for j = 1:m
    b(j) = b(j) + temp;
    temp = temp*xData(i);
  end
  temp = 1;
  for j = 1:2*m-1
    s(j) = s(j) + temp;
    temp = temp*xData(i);
  end
end
for i = 1:m
  for j = 1:m
    A(i,j) = s(i+j-1);
  end
end

% Rearrange coefficients so that coefficient of x^(m-1) is first
coeff = flipdim(gaussPiv(A,b),1);