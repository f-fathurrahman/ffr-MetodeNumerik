func = @(x) ( (1 - x^2)^1.5 )

ExactIntegral = 3*pi/8

for N = 10:10:90
  Integral = gaussQuad( func, -1, 1, N );
  fprintf('%5d  %18.10f %18.10e\n', N, Integral, abs(Integral-ExactIntegral))
end

