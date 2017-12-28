function r = richardson(r,k)
for j = k-1:-1:1
  c = 4^(k-j);
  r(j) = ( c*r(j+1) - r(j) ) / (c-1);
end

