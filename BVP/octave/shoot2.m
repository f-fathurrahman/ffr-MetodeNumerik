
addpath('../../AkarPersamaan/octave')
addpath('../../IVP/octave')

function F = dEqs(x,y)
  % First-order differential
  F = [y(2), -3*y(1)*y(2)]; % equations.
endfunction

function y = inCond(u) % Initial conditions (u is
  y = [0 u]; % the unknown condition).
endfunction


%function shoot2

% Shooting method for 2nd-order boundary value problem in Example 8.1.
xStart = 0;
xStop = 2; % Range of integration.
h = 0.1; % Step size.
freq = 2; % Frequency of printout.
u1 = 1;
u2 = 2; % Trial values of unknown
% initial condition u.
x = xStart;

function r = residual(u)
  % Boundary residual.
  x = 0;
  xStop = 2;
  h = 0.1;
  [xSol,ySol] = runKut4(@dEqs,x,inCond(u),xStop,h);
  r = ySol(size(ySol,1),1) - 1;
endfunction

u = ridder(@residual,u1,u2);
[xSol,ySol] = runKut4(@dEqs,x,inCond(u),xStop,h);
printSol(xSol,ySol,freq)


%endfunction
