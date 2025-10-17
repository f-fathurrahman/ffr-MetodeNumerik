function [Edata] = DOP3Dmudw(dw,input)

dx =input.dx; gam0=input.gamma_0; xR = input.xR; 

delta  = (4*pi/3)^(-1/3) * dx;     % radius sphere with area of dx^3
arg    = gam0*delta; 
factor = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 

Edata = zeros(1,input.NR);

for p = 1 : input.NR     
  % Shifted grid (1) ------------------------------------------------------
    X1  = xR(1,p)-input.X1-dx/2; 
    X2  = xR(2,p)-input.X2; 
    X3  = xR(3,p)-input.X3;
    DIS = sqrt(X1.^2+X2.^2+X3.^2); 
    X1  = X1./DIS;  X2=X2./DIS;  X3= X3./DIS;
      G = factor * exp(-gam0*DIS)./(4*pi*DIS);  dG = -(gam0+1./DIS) .* G;  
   d1_G = X1.*dG;    d2_G = X2.*dG;   d3_G = X3.*dG; 
  E1rfl = -2*dx^2 * sum( d2_G(:).*dw{2,1}(:) + d3_G(:).*dw{3,1}(:));
  E2rfl = -2*dx^2 * sum(-d1_G(:).*dw{2,1}(:)                      );
  E3rfl = -2*dx^2 * sum(                     - d1_G(:).*dw{3,1}(:));
  
  % Shifted grid (2) ------------------------------------------------------
    X1  = xR(1,p)-input.X1; 
    X2  = xR(2,p)-input.X2-dx/2; 
    X3  = xR(3,p)-input.X3; 
    DIS = sqrt(X1.^2+X2.^2+X3.^2);  
    X1  = X1./DIS;  X2=X2./DIS;  X3= X3./DIS;
      G = factor * exp(-gam0*DIS)./(4*pi*DIS);  dG = -(gam0+1./DIS) .* G;  
   d1_G = X1.*dG;    d2_G = X2.*dG;   d3_G = X3.*dG;   
  E1rfl = E1rfl -2*dx^2 * sum(-d2_G(:).*dw{1,2}(:)                      );
  E2rfl = E2rfl -2*dx^2 * sum( d1_G(:).*dw{1,2}(:) + d3_G(:).*dw{3,2}(:));
  E3rfl = E3rfl -2*dx^2 * sum(                     - d2_G(:).*dw{3,2}(:));
  
  % Shifted grid(3) -------------------------------------------------------
     X1 = xR(1,p)-input.X1; 
     X2 = xR(2,p)-input.X2; 
     X3 = xR(3,p)-input.X3-dx/2; 
    DIS = sqrt(X1.^2+X2.^2+X3.^2);  
    X1 = X1./DIS;  X2=X2./DIS;  X3= X3./DIS;
      G = factor * exp(-gam0*DIS)./(4*pi*DIS);  dG = -(gam0+1./DIS) .* G; 
   d1_G = X1.*dG;    d2_G = X2.*dG;   d3_G = X3.*dG;  
  E1rfl = E1rfl -2*dx^2 * sum(-d3_G(:).*dw{1,3}(:)                     );
  E2rfl = E2rfl -2*dx^2 * sum(                     - d3_G(:).*dw{2,3}(:));
  E3rfl = E3rfl -2*dx^2 * sum( d1_G(:).*dw{1,3}(:) + d2_G(:).*dw{2,3}(:));
  
  Edata(1,p) = sqrt(abs(E1rfl)^2  + abs(E2rfl)^2  + abs(E3rfl^2));
  
end % p_loop