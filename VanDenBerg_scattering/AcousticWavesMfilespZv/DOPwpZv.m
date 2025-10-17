function [data] = DOPwpZv(w_p,w_Zv,input)
global nDIM;
gam0 = input.gamma_0;    dx = input.dx;   xR = input.xR;
data = zeros(1,input.NR); 

if nDIM == 1  
  x1     = input.X1; 
  delta  = 0.5 * dx;  factor = sinh(gam0*delta) / (gam0*delta); 
  DIS      = abs(x1-xR(1,1));             % receiver for reflected field
  G        = factor /(2*gam0) * exp(-gam0.*DIS);
  d1_G     = - gam0 * sign(xR(1,1)-x1) .* G;  
 data(1,1) = (gam0^2 * dx) * sum(G .* w_p) ...
               - gam0 * dx * sum(d1_G(:) .* w_Zv{1}(:)); 
  DIS      = abs(x1-xR(1,2));             % receiver for transmitted field
  G        = factor /(2*gam0) * exp(-gam0.*DIS);
  d1_G     = - gam0 * sign(xR(1,2)-x1) .* G; 
 data(1,2) = (gam0^2 * dx) * sum(   G(:) .* w_p(:))     ... 
              - gam0 * dx  * sum(d1_G(:) .* w_Zv{1}(:));
  
elseif nDIM == 2 
  X1    = input.X1;   X2 = input.X2;   
  delta = (pi)^(-1/2)*dx; factor = 2*besseli(1,gam0*delta) / (gam0*delta);  
  for p = 1 : input.NR    
     DIS     = sqrt((X1-xR(1,p)).^2 +(X2-xR(2,p)).^2);
     G       = factor /(2*pi).* besselk(0,gam0*DIS);    
     d_G     = - factor * gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
    d1_G     = ((xR(1,p)-X1)./DIS) .* d_G;
    d2_G     = ((xR(2,p)-X2)./DIS) .* d_G; 
   data(1,p) = (gam0^2 * dx^2) * sum(   G(:) .* w_p(:))     ... 
                 - gam0 * dx^2 * sum(d1_G(:) .* w_Zv{1}(:)) ...
                 - gam0 * dx^2 * sum(d2_G(:) .* w_Zv{2}(:));       
  end % p_loop
 
elseif nDIM == 3
  X1     = input.X1;   X2 = input.X2;   X3 = input.X3;
  delta  = (4*pi/3)^(-1/3) * dx;   arg     = gam0*delta; 
 factor  = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 
  for p = 1 : input.NR
     DIS     = sqrt((X1-xR(1,p)).^2 + (X2-xR(2,p)).^2 + (X3-xR(3,p)).^2);
     G       = factor * exp(-gam0*DIS) ./ (4*pi*DIS); 
     d_G     = (-1./DIS - gam0) .* G;
     d1_G    = ((xR(1,p)-X1)./DIS) .* d_G;
     d2_G    = ((xR(2,p)-X2)./DIS) .* d_G;   
     d3_G    = ((xR(3,p)-X3)./DIS) .* d_G;   
   data(1,p) = (gam0^2 * dx^3) * sum(   G(:) .* w_p(:))     ... 
                 - gam0 * dx^3 * sum(d1_G(:) .* w_Zv{1}(:)) ... 
                 - gam0 * dx^3 * sum(d2_G(:) .* w_Zv{2}(:)) ... 
                 - gam0 * dx^3 * sum(d3_G(:) .* w_Zv{3}(:)); 
  end % p_loop
end % if