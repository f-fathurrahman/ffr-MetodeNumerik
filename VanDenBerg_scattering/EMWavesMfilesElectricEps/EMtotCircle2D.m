clear all; clc; close all; clear workspace
input = initEM();       a = input.a;       gam0 = input.gamma_0; 
eps_sct = input.eps_sct;               mu_sct = input.mu_sct;
gam_sct = gam0 * sqrt(eps_sct*mu_sct); Z_sct  = sqrt(mu_sct/eps_sct);  

% (1) Transform Cartesian coordinates to polar ccordinates  
  xS = input.xS;  rS = sqrt(xS(1)^2+xS(2)^2);    phiS = atan2(xS(2),xS(1)); 
  X1 = input.X1;  X2 = input.X2; 
  R  = sqrt(X1.^2 + X2.^2 + 1e-16);              PHI  = atan2(X2,X1); 

% (2) Compute coefficients of Bessel series expansion ---------------------
  arg0 = gam0 * input.a;  args = gam_sct*input.a; 
  M = 75;                                   % increase M for more accuracy
  A = zeros(1,M); B = zeros(1,M);
  for m = 1 : M                      
    Ib0 = besseli(m,arg0);     dIb0 =  besseli(m+1,arg0) + m/arg0 * Ib0; 
    Ibs = besseli(m,args);     dIbs =  besseli(m+1,args) + m/args * Ibs; 
    Kb0 = besselk(m,arg0);     dKb0 = -besselk(m+1,arg0) + m/arg0 * Kb0; 
    A(m)  = -(Z_sct * dIbs*Ib0 - dIb0*Ibs) / (Z_sct * dIbs*Kb0 - dKb0*Ibs);
    B(m)  =  (Ib0 + A(m) * Kb0) / Ibs;  
  end
  
% (3) Compute exterior and interior  --------------------------------------
    Er_out=zeros(size(R)); Ephi_out=zeros(size(R)); ZH3_out=zeros(size(R));
    Er_int=zeros(size(R)); Ephi_int=zeros(size(R)); ZH3_int=zeros(size(R));
  for m = 1 : M 
    arg0 = gam0*R;      Kb0 = besselk(m,arg0);     
                       dKb0 =-besselk(m+1,arg0) + m./arg0 .* Kb0;
                        Ib0 = besseli(m,arg0);    
                       dIb0 = besseli(m+1,arg0) + m./arg0 .* Ib0;                     
    args = gam_sct*R;   Ibs = besseli(m,args); 
                       dIbs = besseli(m+1,args) + m./args .* Ibs; 
    argS = gam0*rS;     KbS = besselk(m,argS); 
  
    Er_out   = Er_out   + 2*m^2* (Ib0+A(m)*Kb0) .* KbS .* cos(m*(PHI-phiS));
    Ephi_out = Ephi_out - 2*m  *(dIb0+A(m)*dKb0).* KbS .* sin(m*(PHI-phiS));
    ZH3_out  = ZH3_out  - 2*m  * (Ib0+A(m)*Kb0) .* KbS .* sin(m*(PHI-phiS)); 
    Er_int   = Er_int   + 2*m^2*    B(m)* Ibs   .* KbS .* cos(m*(PHI-phiS));
    Ephi_int = Ephi_int - 2*m  *    B(m)*dIbs   .* KbS .* sin(m*(PHI-phiS));
    ZH3_int  = ZH3_int  - 2*m  *    B(m)* Ibs   .* KbS .* sin(m*(PHI-phiS));           
  end % m_loop
 
  Ephi_out(ceil(input.N1/2),ceil(input.N2/2)) = 0;
  ZH3_out(ceil(input.N1/2),ceil(input.N2/2))  = 0;
  Er   =  (                   1/(2*pi) * Er_out  ./R ./rS) .* (R > a)  ...
        + (      1/eps_sct  * 1/(2*pi) * Er_int  ./R ./rS) .* (R < a);
  Ephi =  ( gam0            * 1/(2*pi) * Ephi_out    ./rS) .* (R > a)  ...
        + ( gam_sct/eps_sct * 1/(2*pi) * Ephi_int    ./rS) .* (R < a); 
  ZH3  =  (-gam0            * 1/(2*pi) * ZH3_out     ./rS) .* (R > a)  ...
        + (-gam0            * 1/(2*pi) * ZH3_int     ./rS) .* (R < a);
  
% (4) Transform to Cartesian vectors --------------------------------------
  E{1} = cos(PHI) .* Er  - sin(PHI) .* Ephi;
  E{2} = sin(PHI) .* Er  + cos(PHI) .* Ephi;
  

x1 = X1(:,1); x2 = X2(1,:); phi = 0:.01:2*pi;
   set(figure,'Units','centimeters','Position',[5 5 18 12]);             
subplot(1,3,1)     
IMAGESC(x1,x2,abs(E{1}));  hold on; 
    title('\fontsize{13} abs(E_1)');  plot(a*cos(phi),a*sin(phi),'w'); 
subplot(1,3,2)   
IMAGESC(x1,x2,abs(E{2}));  hold on; 
    title('\fontsize{13} abs(E_2)');  plot(a*cos(phi),a*sin(phi),'w');
subplot(1,3,3)
IMAGESC(x1,x2,abs(ZH3));   hold on; 
title('\fontsize{13} abs(ZH_3)');   plot(a*cos(phi),a*sin(phi),'w'); 


 