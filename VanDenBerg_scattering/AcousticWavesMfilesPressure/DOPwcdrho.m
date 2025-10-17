function [data] = DOPwcdrho(w_c,w_drho,input)
global nDIM;
gam0 = input.gamma_0;    dx = input.dx;   xR = input.xR;
data = zeros(1,input.NR); 

if nDIM == 1   
  x1     = input.X1; 
  delta  = 0.5 * dx;  factor = sinh(gam0*delta) / (gam0*delta); 
  DIS      = abs(x1-xR(1,1));              % receiver for reflected field
  G        = 1/(2*gam0) * exp(-gam0.*DIS);
  d1_G     = - gam0 * sign(xR(1,1)-x1) .* G;  
 data(1,1) = (gam0^2 * dx) * factor * sum(G .* w_c); 
 data(1,1) = data(1,1) - dx * sum(d1_G(:) .* w_drho{1}(:)); 
  DIS      = abs(x1-xR(1,2));            % receiver for transmitted field
  G        = 1/(2*gam0) * exp(-gam0.*DIS);
  d1_G     = - gam0 * sign(xR(1,2)-x1) .* G; 
 data(1,2) = (gam0^2 * dx) * factor * sum(G(:) .* w_c(:)); 
 data(1,2) = data(1,2) - dx * factor * sum(d1_G(:) .* w_drho{1}(:));
  
elseif nDIM == 2 
  X1    = input.X1;   X2 = input.X2;   
  delta = (pi)^(-1/2)*dx; factor = 2*besseli(1,gam0*delta) / (gam0*delta);  
  for p = 1 : input.NR    
     DIS     = sqrt((X1-xR(1,p)).^2 +(X2-xR(2,p)).^2);
     G       = 1/(2*pi).* besselk(0,gam0*DIS);    
     d_G     = - gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
   data(1,p) = (gam0^2 * dx^2) * factor * sum(G(:) .* w_c(:)); 
    d1_G     = ((xR(1,p)-X1)./DIS) .* d_G; 
   data(1,p) = data(1,p) - dx^2 * sum(d1_G(:) .* w_drho{1}(:)); 
    d2_G     = ((xR(2,p)-X2)./DIS) .* d_G; 
   data(1,p) = data(1,p) - dx^2 * sum(d2_G(:) .* w_drho{2}(:));       
  end % p_loop
 
elseif nDIM == 3 
  X1     = input.X1;   X2 = input.X2;   X3 = input.X3;
  delta  = (4*pi/3)^(-1/3) * dx;   arg     = gam0*delta; 
 factor  = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 
  for p = 1 : input.NR
     DIS     = sqrt((X1-xR(1,p)).^2 + (X2-xR(2,p)).^2 + (X3-xR(3,p)).^2);
     G       = exp(-gam0*DIS) ./ (4*pi*DIS); 
     d_G     = (-1./DIS - gam0) .* exp(-gam0*DIS) ./ (4*pi*DIS);
   data(1,p) = (gam0^2 * dx^3) * factor * sum(G(:) .* w_c(:)); 
     d1_G    = ((xR(1,p)-X1)./DIS) .* d_G;   
   data(1,p) = data(1,p) - dx^3 * sum(d1_G(:) .* w_drho{1}(:)); 
     d2_G    = ((xR(2,p)-X2)./DIS) .* d_G; 
   data(1,p) = data(1,p) - dx^3 * sum(d2_G(:) .* w_drho{2}(:));  
     d3_G    = ((xR(3,p)-X3)./DIS) .* d_G; 
   data(1,p) = data(1,p) - dx^3 * sum(d3_G(:) .* w_drho{3}(:)); 
  end % p_loop
end % if
