% Example 2.14 (Matrix inversion)
n = 6;
d = ones(n,1)*2;
e = -ones(n-1,1);
c = e;
d(n) = 5;
[c,d,e] = LUdec3(c,d,e);
for i = 1:n
  b = zeros(n,1);
  b(i) = 1;
  Ainv(:,i) = LUsol3(c,d,e,b);
end
Ainv
