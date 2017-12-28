function I = gaussQuad( func, a, b, N )
%
% func: handle of function to be integrated
% a, b: integration limits
% N: order of integration

% mapping
c1 = ( b + a )/2;
c2 = ( b - a )/2;

% nodes and weights
[x,w] = gaussNodes(N);

s = 0;
for i = 1:length(x)
  y = feval( func, c1 + c2*x(i) ); % value of func at node i
  s = s + w(i)*y;
end
I = c2*s;



