function [data] = DOPdw(dw,data,input)
global nDIM;
gam0 = input.gamma_0; xR   = input.xR;  dx = input.dx;   
   
if nDIM == 1
    
   x1      = input.X1; 
   DIS     = abs(xR(1,1)-(x1+dx/2));       % receiver for reflected field
   G       = 1/(2*gam0) * exp(-gam0.*DIS);
   d1_G    = - gam0 * sign(xR(1,1)-(x1+dx/2)) .* G; 
 data(1,1) = data(1,1) + 2 * sum(d1_G(:) .* dw{1}(:)); 
   DIS     = abs(xR(1,2)-(x1+dx/2));       % receiver for transmitted field
   G       = 1/(2*gam0) * exp(-gam0.*DIS);
   d1_G    = - gam0 * sign(xR(1,2)-(x1+dx/2)) .* G; 
 data(1,2) = data(1,2) + 2 * sum(d1_G(:) .* dw{1}(:));
  
elseif nDIM == 2
    
  X1 = input.X1;   X2 = input.X2;   
  for p = 1 : input.NR    
    DIS      = sqrt((xR(1,p)-(X1+dx/2)).^2 + (xR(2,p)-X2).^2);
    d_G      = - gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
    d1_G     = ((xR(1,p)-(X1+dx/2))./DIS) .* d_G; 
   data(1,p) = data(1,p) + 2 * dx * sum(d1_G(:) .* dw{1}(:)); 
    DIS      = sqrt((xR(1,p)-X1).^2 + (xR(2,p)-(X2+dx/2)).^2);
    d_G      = - gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
    d2_G     = ((xR(2,p)-(X2+dx/2))./DIS) .* d_G; 
   data(1,p) = data(1,p) + 2 * dx * sum(d2_G(:) .* dw{2}(:));       
  end % p_loop
 
elseif nDIM == 3 
    
  X1 = input.X1;   X2 = input.X2;   X3 = input.X3;
 for p = 1 : input.NR    
    DIS     = sqrt((xR(1,p)-(X1+dx/2)).^2+(xR(2,p)-X2).^2+(xR(3,p)-X3).^2);
    d_G     = (-1./DIS-gam0) .* exp(-gam0*DIS) ./ (4*pi*DIS);
    d1_G    = ((xR(1,p)-(X1+dx/2))./DIS) .* d_G;   
  data(1,p) = data(1,p) + 2 * dx^2 * sum(d1_G(:) .* dw{1}(:)); 
    DIS     = sqrt((xR(1,p)-X1).^2+(xR(2,p)-(X2+dx/2)).^2+(xR(3,p)-X3).^2);
    d_G     = (-1./DIS-gam0) .* exp(-gam0*DIS) ./ (4*pi*DIS);
    d2_G    = ((xR(2,p)-(X2+dx/2))./DIS) .* d_G; 
  data(1,p) = data(1,p) + 2 * dx^2 * sum(d2_G(:) .* dw{2}(:));  
    DIS     = sqrt((xR(1,p)-X1).^2+(xR(2,p)-X2).^2+(xR(3,p)-(X3+dx/2)).^2);
    d_G     = (-1./DIS-gam0) .* exp(-gam0*DIS) ./ (4*pi*DIS);
    d3_G    = ((xR(3,p)-(X3+dx/2))./DIS) .* d_G; 
  data(1,p) = data(1,p) + 2 * dx^2 * sum(d3_G(:) .* dw{3}(:)); 
  end % p_loop
  
end % if
