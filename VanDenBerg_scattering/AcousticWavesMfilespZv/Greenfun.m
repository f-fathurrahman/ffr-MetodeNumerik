function [IntG,IntdG] = Greenfun(temp,Sign,input)
global nDIM; 

gam0 = input.gamma_0;
dx   = input.dx; 
if nDIM == 1     
   delta   = 0.5 * dx; 
   factor  = sinh(gam0*delta) / (gam0*delta); 
   X1      = temp.X1fft;  
   DIS     = abs(X1);
   G       = 1/(2*gam0) * exp(-gam0.*DIS) * factor^2;
   IntG    = (gam0^2 * dx) * G;                 % integral includes gam0^2
 IntG(1)   = 1 - exp(-gam0*delta) * factor;  %----------------------------
   d_G     = - gam0 * G;
 d1_G      = Sign.X1fft .* d_G;
 IntdG{1}  = dx * d1_G;
   
elseif nDIM == 2   
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
  
elseif nDIM == 3     % integral over sphere with equivalent volume of dx^3   
   delta   = (4*pi/3)^(-1/3) * dx;       % radius sphere with area of dx^3
   arg     = gam0*delta;
   factor  = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2;              
    X1     = temp.X1fft;  X2 = temp.X2fft;  X3 = temp.X3fft;
    DIS    = sqrt(X1.^2 + X2.^2 + X3.^2);
DIS(1,1,1) = 1;                    % avoid Green's singularity for DIS = 0  
   G       = exp(-gam0*DIS) ./ (4*pi*DIS) * factor^2;  
   IntG    = (gam0^2 * dx^3) * G;               % integral includes gam0^2
IntG(1,1,1)= 1 - (1+gam0*delta) * exp(-gam0*delta) * factor;
  d_G      = (-1./DIS - gam0) .* G;
  d1_G     = (Sign.X1fft .* X1./DIS) .* d_G; 
  d2_G     = (Sign.X2fft .* X2./DIS) .* d_G;  
  d3_G     = (Sign.X3fft .* X3./DIS) .* d_G; 
  IntdG{1} = dx^3 * d1_G;    IntdG{1}(1,1,1) = 0;
  IntdG{2} = dx^3 * d2_G;    IntdG{2}(1,1,1) = 0;
  IntdG{3} = dx^3 * d3_G;    IntdG{3}(1,1,1) = 0;
end