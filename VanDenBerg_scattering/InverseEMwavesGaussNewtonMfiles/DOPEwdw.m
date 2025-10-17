function [Ercv] = DOPEwdw(w_E,dwE,input)

gam0 = input.gamma_0;      dx = input.dx;   xR = input.xR; 
 
delta  = (pi)^(-1/2)*dx;            % radius circle with area dx^2   
factor = 2 * besseli(1,gam0*delta) / (gam0*delta); 
 
for p = 1 : input.NR   
   % Non shifted grid(0)  
       X1  = xR(1,p)-input.X1;     X2 = xR(2,p)-input.X2; 
     DIS  = sqrt(X1.^2 + X2.^2);   X1 = X1./DIS;  X2 = X2./DIS;
      G   =   factor * 1/(2*pi).* besselk(0,gam0*DIS); 
     dG   = - factor * gam0 * 1/(2*pi).* besselk(1,gam0*DIS);  
    d1_G  =  X1 .* dG;      
    d2_G  =  X2 .* dG;  
    E1rfl = gam0^2 * dx^2 * sum( G(:) .* w_E{1}(:) ); 
    E2rfl = gam0^2 * dx^2 * sum( G(:) .* w_E{2}(:) ); 

   % Shifted grid(1)
      X1  = xR(1,p) - input.X1 - dx/2;            X2 = xR(2,p)-input.X2; 
     DIS  = sqrt(X1.^2 + X2.^2);   X1 = X1./DIS;  X2 = X2./DIS;
     dG   = - factor * gam0 * 1/(2*pi).* besselk(1,gam0*DIS);
    d1_G  =  X1 .* dG; 
    d2_G  =  X2 .* dG; 
    E1rfl = E1rfl + 2 * dx * sum(d1_G(:) .* dwE{1}(:));
    E2rfl = E2rfl + 2 * dx * sum(d2_G(:) .* dwE{1}(:));
    
   % Shifted grid(2)
      X1  = xR(1,p) - input.X1;    X2 = xR(2,p) - input.X2 - dx/2; 
     DIS  = sqrt(X1.^2 + X2.^2);   X1 = X1./DIS;  X2 = X2./DIS;
     dG   = - factor * gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);   
    d1_G  =  X1 .* dG; 
    d2_G  =  X2 .* dG;  
   E1rfl  = E1rfl + 2 * dx * sum(d1_G(:) .* dwE{2}(:));  
   E2rfl  = E2rfl + 2 * dx * sum(d2_G(:) .* dwE{2}(:));  

   Ercv{1}(p,1) = E1rfl;
   Ercv{2}(p,1) = E2rfl;
 end % p_loop