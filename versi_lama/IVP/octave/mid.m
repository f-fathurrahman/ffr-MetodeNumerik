function y = mid(dEqs,x,y,xStop,nSteps)
  % Midpoint formulas.
  h = (xStop - x)/nSteps;
  y0 = y;
  y1 = y0 + h*feval(dEqs,x,y0);
  for i = 1:nSteps-1
    x = x + h;
    y2 = y0 + 2.0*h*feval(dEqs,x,y1);
    y0 = y1;
    y1 = y2;
  end
  y = 0.5*(y1 + y0 + h*feval(dEqs,x,y2));
endfunction