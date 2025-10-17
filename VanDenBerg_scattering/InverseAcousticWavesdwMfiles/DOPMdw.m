function [data] = DOPMdw(dw,input)
gam0 = input.gamma_0;    dx = input.dx;   xR = input.xR;
X1   = input.X1;         X2 = input.X2;   

data = zeros(input.NR,1); 
for p = 1 : input.NR    
   DIS      = sqrt((xR(1,p)-(X1+dx/2)).^2 + (xR(2,p)-X2).^2);
   d_G      = - gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
   d1_G     = ((xR(1,p)-(X1+dx/2))./DIS) .* d_G; 
  data(p,1) = data(p,1) + 2 * dx * sum(d1_G(:) .* dw{1}(:)); 
   
   DIS      = sqrt((xR(1,p)-X1).^2 + (xR(2,p)-(X2+dx/2)).^2);
   d_G      = - gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
   d2_G     = ((xR(2,p)-(X2+dx/2))./DIS) .* d_G; 
  data(p,1) = data(p,1) + 2 * dx * sum(d2_G(:) .* dw{2}(:));  
end % p_loop