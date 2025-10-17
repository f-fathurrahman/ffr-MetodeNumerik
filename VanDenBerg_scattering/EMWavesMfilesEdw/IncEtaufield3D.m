function [Einc_tau] = IncEtaufield3D(input) 

gam0 = input.gamma_0;  xS = input.xS;    dx = input.dx;
delta  = (4*pi/3)^(-1/3) * dx;     % radius sphere with area of dx^3
arg    = gam0*delta; factor = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 
  
% Shifted grid (1) --------------------------------------------------------
   X1   = input.X1+dx/2-xS(1);  
   X2   = input.X2-xS(2);  
   X3   = input.X3-xS(3); 
   DIS  = sqrt(X1.^2 + X2.^2 + X3.^2); 
   X1   = X1./DIS;             
   X2   = X2./DIS;        
   X3   = X3./DIS;
    G   = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  
   dG21 = ( 3*X2.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2*X2.*X1 ) .* G;
   dG31 = ( 3*X3.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2*X3.*X1 ) .* G;
   
   Einc_tau{2,1} = - dG21; 
   Einc_tau{3,1} = - dG31;
   
 % Shifted grid (2) -------------------------------------------------------
   X1   = input.X1-xS(1); 
   X2   = input.X2+dx/2-xS(2); 
   X3   = input.X3-xS(3); 
   DIS  = sqrt(X1.^2 + X2.^2 + X3.^2); 
   X1   = X1./DIS;   
   X2   = X2./DIS;             
   X3   = X3./DIS;
     G  = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  
   dG11 = ((3*X1.*X1 -1).* (gam0./DIS + 1./DIS.^2) + gam0^2*X1.*X1 ) .* G;
   dG31 = ( 3*X3.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2*X3.*X1 ) .* G;

   Einc_tau{1,2} =  - (-gam0^2 * G + dG11);  
   Einc_tau{3,2} =  - dG31;
   
 % Shifted grid (3) -------------------------------------------------------
  X1    = input.X1-xS(1);  
  X2    = input.X2-xS(2);  
  X3    = input.X3+dx/2-xS(3); 
  DIS   = sqrt(X1.^2 + X2.^2 + X3.^2); 
  X1    = X1./DIS;         
  X2    = X2./DIS;         
  X3    = X3./DIS;
    G   = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  
   dG11 = ((3*X1.*X1 -1).* (gam0./DIS + 1./DIS.^2) + gam0^2*X1.*X1 ) .* G;
   dG21 = ( 3*X2.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2*X2.*X1 ) .* G;
 
   Einc_tau{1,3} = - (-gam0^2 * G + dG11);  
   Einc_tau{2,3} = - dG21; 