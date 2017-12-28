% Example 6.11 (Gauss-Legendre quadrature)
func = @(x) ((sin(x)/x)^2);

a = 0;
b = pi;
Iexact = 1.41815;

for n = 2:12
  I = gaussQuad(func, a, b, n);
  diffI = abs(I - Iexact);
  fprintf('%5d  %18.10f %18.10e\n', n, I, diffI);
  if diffI < 1e-6
    fprintf('Convergence achieved for n = %d\n', n)
    fprintf('I = %18.10f\n', I)
    break
  end
end

