function Av = fex2_18(v)
% Computes the product A*v in Example 2.18
n = length(v);
Av = zeros(n,1);
Av(1) = 2*v(1) - v(2) + v(n);
Av(2:n-1) = -v(1:n-2) + 2*v(2:n-1) - v(3:n);
Av(n) = -v(n-1) + 2*v(n) + v(1);