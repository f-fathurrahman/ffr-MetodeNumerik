// evaluate polynomial
// p = a(1)*x^n + a(2)*x^(n-1) + ... + a(n+1)
// and its first two derivatives dp and d_dp

function [p, dp, d_dp] = evalpoly(a,x)

  n = length(a) - 1

  p = a(1)
  dp = 0.0
  d_dp = 0.0

  for i = 1:n
    d_dp = d_dp*x + 2.0*dp
    dp = dp*x + p
    p = p*x + a(i+1)
  end

endfunction