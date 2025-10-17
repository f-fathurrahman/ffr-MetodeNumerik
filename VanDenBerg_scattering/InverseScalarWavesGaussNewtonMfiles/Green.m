function [IntG] = Green(temp,input)

  dx = input.dx; 
gam0 = input.gamma_0;

 X1      = temp.X1fft;  
 X2      = temp.X2fft;  
DIS      = sqrt(X1.^2 + X2.^2); 
DIS(1,1) = 1;                     % avoid Green's singularity for DIS = 0
 G       = 1/(2*pi).* besselk(0,gam0*DIS);
 delta   = (pi)^(-1/2) * dx;      % radius circle with area of dx^2
 factor  = 2 * besseli(1,gam0*delta) / (gam0*delta);   
  
% Integral of Green function including gam0^2
IntG      = (gam0^2 * dx^2) * factor^2 * G;   
IntG(1,1) = 1 - gam0*delta * besselk(1,gam0*delta) * factor; 