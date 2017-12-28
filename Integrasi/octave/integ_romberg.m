function [I, numEval] = integ_romberg( func, a, b, tol, kMax )
% 
% func: handle of function to be integrated
% a, b: limits of integration
% tol: error tolerance
% kMax: limit on the number of panel doublings

%
% set up default arguments
%
if nargin < 5
  kMax = 20;
end

if nargin < 4
  tol = 1.e4 * eps;
end

r    = zeros(kMax);
r(1) = integ_trapezoid_recursive( func, a, b, 0, 1 );
rOld = r(1);

for k = 2:kMax
  r(k) = integ_trapezoid_recursive( func, a, b, r(k-1), k );
  r = richardson(r,k);
  if abs( r(1) - rOld ) < tol
    numEval = 2^(k-1) + 1;
    I = r(1);
    return
  end
  rOld = r(1);
end
error('integ_romberg: ERROR: Failed to converge')

