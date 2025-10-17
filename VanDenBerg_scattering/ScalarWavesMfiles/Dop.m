function [data] = Dop(w,input)
global nDIM;
gam0 = input.gamma_0;    
dx   = input.dx;   
xR   = input.xR;

data = zeros(1,input.NR); 

if nDIM == 1  

  x1      = input.X1; 
  delta   = 0.5 * dx;
  factor  = sinh(gam0*delta) / (gam0*delta); 
  DIS     = abs(xR(1,1)-x1);              % receiver for reflected field
  G       = 1/(2*gam0) * exp(-gam0.*DIS); 
 data(1,1)= (gam0^2 * dx) * factor * sum(G.*w); 
  DIS     = abs(xR(1,2)-x1);              % receiver for transmitted field
  G       = 1/(2*gam0) * exp(-gam0.*DIS);
 data(1,2)= (gam0^2 * dx) * factor * sum(G.*w); 
  
elseif nDIM == 2 

  X1      = input.X1; X2 = input.X2;   
  delta   = (pi)^(-1/2) * dx;             % radius circle with area of dx^2 
  factor  = 2 * besseli(1,gam0*delta) / (gam0*delta);  
  for p = 1 : input.NR    
     DIS     = sqrt((xR(1,p)-X1).^2 +(xR(2,p)-X2).^2);
     G       = 1/(2*pi).* besselk(0,gam0*DIS);
   data(1,p) = (gam0^2 * dx^2) * factor * sum(G(:).*w(:));       
  end % p_loop
 
elseif nDIM == 3

  X1      = input.X1; X2 = input.X2; X3 = input.X3;
  delta   = (4*pi/3)^(-1/3) * dx;         % radius sphere with area of dx^3
  arg     = gam0*delta;
  factor  = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 
  for p = 1 : input.NR
     DIS     = sqrt((xR(1,p)-X1).^2 + (xR(2,p)-X2).^2 + (xR(3,p)-X3).^2);
     G       = exp(-gam0*DIS) ./ (4*pi*DIS);
   data(1,p) = (gam0^2 * dx^3) * factor * sum(G(:).*w(:));      
  end % p_loop
  
end % if