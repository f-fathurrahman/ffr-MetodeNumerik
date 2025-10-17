function [E_inc,dEinc] = IncEfield3D(input) 
gam0   = input.gamma_0;  xS = input.xS;    dx = input.dx;
delta  = (4*pi/3)^(-1/3) * dx;     % radius sphere with area of dx^3
arg    = gam0*delta; factor = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 
  
% Non shifted grid (0) ----------------------------------------------------
   X1   = input.X1-xS(1);   
   X2   = input.X2-xS(2);     
   X3   = input.X3-xS(3); 
   DIS  = sqrt(X1.^2 + X2.^2 + X3.^2); 
   X1   = X1./DIS;          
   X2   = X2./DIS;            
   X3   = X3./DIS;
    G   = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  
   dG11 = ((3*X1.*X1 -1).* (gam0./DIS + 1./DIS.^2) + gam0^2*X1.*X1 ) .* G;
   dG21 = ( 3*X2.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2*X2.*X1 ) .* G;
   dG31 = ( 3*X3.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2*X3.*X1 ) .* G;
E_inc{1} = - (-gam0^2 * G + dG11);
E_inc{2} = - dG21; 
E_inc{3} = - dG31;
 
% Shifted grid (1) --------------------------------------------------------
   X1   = input.X1+dx/2-xS(1); 
   X2   = input.X2-xS(2);  
   X3   = input.X3-xS(3); 
   DIS  = sqrt(X1.^2 + X2.^2 + X3.^2); 
   X1   = X1./DIS;              
    G   = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  
   dG11 = ((3*X1.*X1 -1).* (gam0./DIS + 1./DIS.^2) + gam0^2*X1.*X1 ) .* G;
 dEinc{1} = - (-gam0^2 * G + dG11);
 % Shifted grid(2) --------------------------------------------------------
   X1   = input.X1-xS(1); 
   X2   = input.X2+dx/2-xS(2); 
   X3   = input.X3-xS(3); 
   DIS  = sqrt(X1.^2 + X2.^2 + X3.^2); 
   X1   = X1./DIS;        
   X2   = X2./DIS;           
     G  = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  
   dG21 = ( 3*X2.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2*X2.*X1 ) .* G;
 dEinc{2} = - dG21; 
 % Shifted grid(3) --------------------------------------------------------
  X1    = input.X1-xS(1);  
  X2    = input.X2-xS(2);  
  X3    = input.X3+dx/2-xS(3); 
  DIS   = sqrt(X1.^2 + X2.^2 + X3.^2); 
  X1    = X1./DIS;                  
  X3    = X3./DIS;
    G   = factor * exp(-gam0*DIS) ./ (4*pi*DIS);  
   dG31 = ( 3*X3.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2*X3.*X1 ) .* G;
  dEinc{3} = - dG31;