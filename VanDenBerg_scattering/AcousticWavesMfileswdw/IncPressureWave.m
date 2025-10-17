function [p_inc,Pinc] = IncPressureWave(input) 
global nDIM;
gam0 = input.gamma_0;  
xS   = input.xS;    dx = input.dx;

if nDIM == 1                       % incident wave on one-dimensional grid  

% Weak form
  delta  = 0.5 * dx;
  factor = sinh(gam0*delta) / (gam0*delta); 
  x1     = input.X1; 
  DIS    = abs(x1-xS(1));           
  p_inc  = factor * 1/(2*gam0) * exp(-gam0.*DIS);  % Grid (0) on midpoints
% Strong form
  DIS     = abs(x1+dx/2-xS(1));
  Pinc{1} = 1/(2*gam0) * exp(-gam0.*DIS);     % Grid (1) on shifted points               

elseif nDIM == 2                   % incident wave on two-dimensional grid

% Weak form
  delta  = (pi)^(-1/2) * dx;       % radius circle with area of dx^2 
  factor = 2 * besseli(1,gam0*delta) / (gam0*delta); 
  X1     = input.X1;   X2 = input.X2; 
  DIS    = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 );    
  p_inc  = factor * 1/(2*pi).* besselk(0,gam0*DIS);% Grid (0) on midpoints 
% Strong form
  DIS     = sqrt( (X1+dx/2-xS(1)).^2 + (X2-xS(2)).^2 );
  Pinc{1} = 1/(2*pi).* besselk(0,gam0*DIS);   % Grid (1) on shifted points               
  DIS     = sqrt( (X1-xS(1)).^2 + (X2+dx/2-xS(2)).^2 );
  Pinc{2} = 1/(2*pi).* besselk(0,gam0*DIS);   % Grid (2) on shifted points
  
elseif nDIM == 3                % incident wave on three-dimensional grid

% Weak form
  delta  = (4*pi/3)^(-1/3) * dx; % radius sphere with area of dx^3
  arg    = gam0*delta;
  factor = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 
  X1     = input.X1;   X2 = input.X2;    X3 =input.X3;
  DIS    = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 + (X3-xS(3)).^2 );
  p_inc  = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  % Grid (0) on midpoints 
% Strong form
  DIS     = sqrt( (X1+dx/2-xS(1)).^2 + (X2-xS(2)).^2 + (X3-xS(3)).^2 );
  Pinc{1} = exp(-gam0*DIS) ./ (4*pi*DIS);     % Grid (1) on shifted points                             
  DIS     = sqrt( (X1-xS(1)).^2 + (X2+dx/2-xS(2)).^2 + (X3-xS(3)).^2 );
  Pinc{2} = exp(-gam0*DIS) ./ (4*pi*DIS);     % Grid (2) on shifted points
  DIS     = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 + (X3+dx/2-xS(3)).^2 );
  Pinc{3} = exp(-gam0*DIS) ./ (4*pi*DIS);     % Grid (3) on shifted points 

end % if