function [u_inc] = IncWave(input) 
global nDIM;

gam0 = input.gamma_0;  
xS   = input.xS;  
dx   = input.dx;

if nDIM == 1                     % incident wave on one-dimensional grid  
 
  x1      = input.X1; 
  DIS     = abs(x1-xS(1));
  G       = 1/(2*gam0) * exp(-gam0.*DIS); 
  delta   = 0.5 * dx;
  factor  = sinh(gam0*delta) / (gam0*delta);                    
  
elseif nDIM == 2                 % incident wave on two-dimensional grid

   X1     = input.X1; X2 = input.X2; 
   DIS    = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 );
   G      = 1/(2*pi).* besselk(0,gam0*DIS);
   delta  = (pi)^(-1/2) * dx;    % radius circle with area of dx^2 
   factor = 2 * besseli(1,gam0*delta) / (gam0*delta);               

elseif nDIM == 3                 % incident wave on three-dimensional grid

   X1     = input.X1; X2 = input.X2; X3 =input.X3;
   DIS    = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 + (X3-xS(3)).^2 );
   G      = exp(-gam0*DIS) ./ (4*pi*DIS);  
   delta  = (4*pi/3)^(-1/3) * dx; % radius sphere with area of dx^3
   arg    = gam0*delta;
   factor = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2;      
           
end

u_inc = factor * G;               % factor for weak form if DIS > delta
