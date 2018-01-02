% Example 2.11 (Solution of pentadiagonal eqs.)
n = 10;
d = 6*ones(n,1);
d(n) = 7;
e = -4*ones(n-1,1);
f = ones(n-2,1);
b = zeros(n,1);
b(1) = 3;
b(n) = 4;
[d,e,f] = LUdec5(d,e,f);
x = LUsol5(d,e,f,b)