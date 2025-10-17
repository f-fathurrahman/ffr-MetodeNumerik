function [IntG] = Green(temp,input)
global nDIM; dx = input.dx; gam0 = input.gamma_0;

if nDIM == 1  
    
   X1      = temp.X1fft;  
   DIS     = abs(X1);
   G       = 1/(2*gam0) * exp(-gam0.*DIS);
   delta   = 0.5 * dx; 
   factor  = sinh(gam0*delta) / (gam0*delta);    
   IntG    = (gam0^2 * dx) * factor^2 * G;   % integral includes gam0^2
   IntG(1) = 1 - exp(-gam0*delta) * factor;

elseif nDIM == 2
    
   X1      = temp.X1fft;  X2 = temp.X2fft;  
   DIS     = sqrt(X1.^2 + X2.^2); 
 DIS(1,1)  = 1;                     % avoid Green's singularity for DIS = 0
   G       = 1/(2*pi).* besselk(0,gam0*DIS);
   delta   = (pi)^(-1/2) * dx;            % radius circle with area of dx^2
   factor  = 2 * besseli(1,gam0*delta) / (gam0*delta);   
   IntG    = (gam0^2 * dx^2) * factor^2 * G; % integral includes gam0^2
 IntG(1,1) = 1 - gam0*delta * besselk(1,gam0*delta) * factor; 
 
elseif nDIM == 3    
    
    X1     = temp.X1fft;  X2 = temp.X2fft;  X3 = temp.X3fft;
    DIS    = sqrt(X1.^2 + X2.^2 + X3.^2);
DIS(1,1,1) = 1;                     % avoid Green's singularity for DIS = 0  
   G       = exp(-gam0*DIS) ./ (4*pi*DIS);  
   delta   = (4*pi/3)^(-1/3) * dx;        % radius sphere with area of dx^3
   arg     = gam0*delta;
   factor  = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2;              
   IntG    = (gam0^2 * dx^3) * factor^2 * G;  % integral includes gam0^2
IntG(1,1,1)= 1 - (1+gam0*delta)*exp(-gam0*delta) * factor;

end