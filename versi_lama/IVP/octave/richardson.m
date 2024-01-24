function r = richardson(r,k,n)
  % Richardson extrapolation.
  for j = k-1:-1:1
    c =(k/(k-1))^(2*(k-j));
    r(j,1:n) =(c*r(j+1,1:n) - r(j,1:n))/(c - 1.0);
  end
endfunction