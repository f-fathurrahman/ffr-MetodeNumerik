function [IntG,IntdG] = Greenfun(temp,Sign,input)

gam0 = input.gamma_0;
dx   = input.dx; 

  delta    = (pi)^(-1/2) * dx;          % radius circle with area of dx^2
  factor   = 2 * besseli(1,gam0*delta) / (gam0*delta);  
  X1       = temp.X1fft;  X2 = temp.X2fft;  
  DIS      = sqrt(X1.^2 + X2.^2); 
DIS(1,1)   = 1;                   % avoid Green's singularity for DIS = 0
   G       = 1/(2*pi).* besselk(0,gam0*DIS) * factor^2;
   IntG    = (gam0^2 * dx^2) * G;              % integral includes gam0^2
IntG(1,1)  = 1 - gam0*delta * besselk(1,gam0*delta) * factor; %---------- 
  d_G      = - gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS) * factor^2; 
  d1_G     = (Sign.X1fft .* X1./DIS) .* d_G;  
  d2_G     = (Sign.X2fft .* X2./DIS) .* d_G;      
IntdG{1}   = dx^2 * d1_G;    IntdG{1}(1,1) = 0;  
IntdG{2}   = dx^2 * d2_G;    IntdG{2}(1,1) = 0;
