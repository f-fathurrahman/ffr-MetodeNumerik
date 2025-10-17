function [p_inc,Zv_inc] = IncAcousticWave(input) 
global nDIM;
gam0 = input.gamma_0;  
xS   = input.xS;    dx = input.dx;

if nDIM == 1                      % incident wave on one-dimensional grid  

% Weak form
  delta     = 0.5 * dx;
  factor    = sinh(gam0*delta) / (gam0*delta); 
  x1        = input.X1; 
  DIS       = abs(x1-xS(1));           
  p_inc     = factor * 1/(2*gam0) * exp(-gam0.*DIS);
  d_p_inc   = - gam0 * p_inc;
  Zv_inc{1} = - (1/gam0) * sign(x1-xS(1)) .* d_p_inc;

elseif nDIM == 2                  % incident wave on two-dimensional grid

% Weak form
  delta     = (pi)^(-1/2) * dx;       % radius circle with area of dx^2 
  factor    = 2 * besseli(1,gam0*delta) / (gam0*delta); 
  X1        = input.X1;   X2 = input.X2; 
  DIS       = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 );    
  p_inc     = factor * 1/(2*pi).* besselk(0,gam0*DIS);
  d_p_inc   = - gam0 * factor * 1/(2*pi).* besselk(1,gam0*DIS);
  Zv_inc{1} = - (1/gam0) * (X1-xS(1))./DIS .* d_p_inc;
  Zv_inc{2} = - (1/gam0) * (X2-xS(2))./DIS .* d_p_inc;
  
elseif nDIM == 3                % incident wave on three-dimensional grid

% Weak form
  delta     = (4*pi/3)^(-1/3) * dx; % radius sphere with area of dx^3
  arg       = gam0*delta;
  factor    = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 
  X1        = input.X1;   X2 = input.X2;    X3 =input.X3;
  DIS       = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 + (X3-xS(3)).^2 );
  p_inc     = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  
  d_p_inc   = - (1./DIS + gam0) .* p_inc;
  Zv_inc{1} = - (1/gam0) * (X1-xS(1))./DIS .* d_p_inc;
  Zv_inc{2} = - (1/gam0) * (X2-xS(2))./DIS .* d_p_inc;
  Zv_inc{3} = - (1/gam0) * (X3-xS(3))./DIS .* d_p_inc;
  
end % if