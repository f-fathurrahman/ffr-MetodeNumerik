function k = splineCurv(xData,yData)
% Returns curvatures of a cubic spline at the knots.
% USAGE: k = splineCurv(xData,yData)
% xData = x-coordinates of data points.
% yData = y-coordinates of data points.
n = length(xData);
c = zeros(n-1,1);
d = ones(n,1);
e = zeros(n-1,1);
k = zeros(n,1);

c(1:n-2) = xData(1:n-2) - xData(2:n-1);
d(2:n-1) = 2*(xData(1:n-2) - xData(3:n));
e(2:n-1) = xData(2:n-1) - xData(3:n);
k(2:n-1) =  6*(yData(1:n-2) - yData(2:n-1))./(xData(1:n-2) - xData(2:n-1)) ... 
          - 6*(yData(2:n-1) - yData(3:n)) ./(xData(2:n-1) - xData(3:n));
[c,d,e] = LUdec3(c,d,e);
k = LUsol3(c,d,e,k);