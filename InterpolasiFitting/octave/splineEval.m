function y = splineEval(xData,yData,k,x)
% Returns value of cubic spline interpolant at x.
% USAGE: y = splineEval(xData,yData,k,x)
% xData = x-coordinates of data points.
% yData = y-coordinates of data points.
% k = curvatures of spline at the knots;
%     returned by the function splineCurv.

i = findSeg(xData,x);

h = xData(i) - xData(i+1);
y =   ((x - xData(i+1))^3/h - (x - xData(i+1))*h)*k(i)/6.0 ...
    - ((x - xData(i))^3/h - (x - xData(i))*h)*k(i+1)/6.0 ...
    + yData(i)*(x - xData(i+1))/h - yData(i+1)*(x - xData(i))/h;

endfunction

function i = findSeg(xData,x)
% Returns index of segment containing x.
iLeft = 1;
iRight = length(xData);
while 1
  if(iRight - iLeft) <= 1
    i = iLeft;
    return
  end
  i = fix((iLeft + iRight)/2);
  if x < xData(i)
    iRight = i;
  else
    iLeft = i;
  end
end

endfunction