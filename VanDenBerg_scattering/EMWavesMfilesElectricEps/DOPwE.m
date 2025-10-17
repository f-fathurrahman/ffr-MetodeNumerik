function [Edata,Hdata] = DOPwE(w_E,input)
global nDIM;
gam0 = input.gamma_0;    dx = input.dx;   xR = input.xR; 

Edata = zeros(1,input.NR); 
Hdata = zeros(1,input.NR); 

if nDIM == 2  
 % Weak form   
 delta  = (pi)^(-1/2)*dx;             % radius circle with area of dx^2
 factor = 2*besseli(1,gam0*delta) / (gam0*delta); 
 
 for p = 1 : input.NR     
      X1  = xR(1,p)-input.X1;     X2 = xR(2,p)-input.X2; 
      DIS = sqrt(X1.^2 + X2.^2);  X1 = X1./DIS; X2 = X2./DIS;
      G   =    factor        * 1/(2*pi).* besselk(0,gam0*DIS); 
 
     dG   =  - factor * gam0 * 1/(2*pi).* besselk(1,gam0*DIS);  
    d1_G  =  X1 .* dG;     d2_G  =  X2 .* dG; 
     
     dG11 = (2*X1.*X1 - 1) .* (-dG./DIS) + gam0^2*X1.*X1 .* G;     
     dG22 = (2*X2.*X2 - 1) .* (-dG./DIS) + gam0^2*X2.*X2 .* G;
     dG21 =  2*X2.*X1      .* (-dG./DIS) + gam0^2*X2.*X1 .* G;
 
    E1rfl = dx^2 * sum( (gam0^2*G(:) - dG11(:)) .* w_E{1}(:)  ...
                                     - dG21(:)  .* w_E{2}(:) ); 
    E2rfl = dx^2 * sum(              - dG21(:)  .* w_E{1}(:)  ...
                       +(gam0^2*G(:) - dG22(:)) .* w_E{2}(:) ); 
   
   ZH3rfl = gam0 * dx^2 * sum( d2_G(:).*w_E{1}(:) - d1_G(:).*w_E{2}(:) );                 
  
   Edata(1,p) = sqrt(abs(E1rfl)^2 + abs(E2rfl)^2);
   Hdata(1,p) = abs(ZH3rfl);   
 
 end; % p_loop
 
 
elseif nDIM == 3 
    
 [Edata,Hdata] = DOP3DwE(w_E,input);
    
end % if