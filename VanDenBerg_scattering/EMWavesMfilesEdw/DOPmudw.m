function [Edata] = DOPmudw(dw,input)
global nDIM;

dx =input.dx; gam0=input.gamma_0; xR = input.xR; 

Edata = zeros(1,input.NR);

if nDIM == 2   
  delta  = (pi)^(-1/2) * dx;         % radius circle with area of dx^2 
  factor = 2 * besseli(1,gam0*delta) / (gam0*delta); 
  
   for p = 1 : input.NR  
   % Shifted grid (1) -----------------------------------------------------
     X1  = xR(1,p)-input.X1-dx/2;    
     X2  = xR(2,p)-input.X2;
     DIS = sqrt(X1.^2 + X2.^2);    
     X1  = X1./DIS;     X2 = X2./DIS;
     d_G = - factor * gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
    d1_G = X1 .* d_G;  d2_G  = X2 .* d_G; 
    
   E1rfl = - 2 * dx * sum(d2_G(:) .* dw{2,1}(:));
   E2rfl =   2 * dx * sum(d1_G(:) .* dw{2,1}(:));
    
   % Shifted grid (2) -----------------------------------------------------
     X1  = xR(1,p)-input.X1;        
     X2  = xR(2,p)-input.X2-dx/2;
     DIS = sqrt(X1.^2 + X2.^2);      
     X1  = X1./DIS;     X2 = X2./DIS;
     d_G = - factor * gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
    d1_G = X1 .* d_G;  d2_G  = X2 .* d_G; 
    
   E1rfl = E1rfl + 2 * dx * sum(d2_G(:) .* dw{1,2}(:));
   E2rfl = E2rfl - 2 * dx * sum(d1_G(:) .* dw{1,2}(:));  
  
   Edata(1,p) = sqrt(abs(E1rfl)^2  + abs(E2rfl)^2);
   end % p_loop
  
elseif nDIM == 3  
 [Edata] = DOP3Dmudw(dw,input);
end % if