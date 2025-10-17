function [data] = HDOPwE(w_E,input)
global nDIM;
gam0 = input.gamma_0;    dx = input.dx;   xR = input.xR; 

data = zeros(1,input.NR); 

if nDIM == 2; 
 delta  = (pi)^(-1/2)*dx;               % radius circle with area dx^2
 factor = 2*besseli(1,gam0*delta) / (gam0*delta); 
 for p = 1 : input.NR     
      X1  = xR(1,p)-input.X1;     X2 = xR(2,p)-input.X2; 
      DIS = sqrt(X1.^2 + X2.^2);  X1 = X1./DIS; X2 = X2./DIS;
      G   =    factor        * 1/(2*pi).* besselk(0,gam0*DIS); 
     dG   =  - factor * gam0 * 1/(2*pi).* besselk(1,gam0*DIS);  
    d1_G  =  X1 .* dG; 
    d2_G  =  X2 .* dG;  
   ZH3rfl = gam0 * dx^2 * sum( d2_G(:).*w_E{1}(:) - d1_G(:).*w_E{2}(:) );
   data(1,p) =  ZH3rfl;
 end % p_loop
  
elseif nDIM == 3; 
 delta  = (4*pi/3)^(-1/3) * dx;   arg     = gam0*delta; 
 factor  = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 
 for p = 1 : input.NR
    X1  = xR(1,p)-input.X1;  X2 = xR(2,p)-input.X2; X3 = xR(3,p)-input.X3;
    DIS = sqrt(X1.^2+X2.^2+X3.^2);  
    X1  = X1./DIS;  X2 = X2./DIS;  X3 = X3./DIS;
    G   = factor * exp(-gam0*DIS) ./ (4*pi*DIS); 
    dG  =  -  (gam0 + 1./DIS) .* G;  
   d1_G = X1 .* dG; 
   d2_G = X2 .* dG; 
   d3_G = X3 .* dG;
   ZH1rfl = -gam0 * dx^3 * sum( d2_G(:).*w_E{3}(:) - d3_G(:).*w_E{2}(:) );
   ZH2rfl = -gam0 * dx^3 * sum( d3_G(:).*w_E{1}(:) - d1_G(:).*w_E{3}(:) );
   ZH3rfl = -gam0 * dx^3 * sum( d1_G(:).*w_E{2}(:) - d2_G(:).*w_E{1}(:) );
   data(1,p) =  sqrt(abs(ZH1rfl)^2 + abs(ZH2rfl)^2 + abs(ZH3rfl^2));
 end % p_loop
end % if
