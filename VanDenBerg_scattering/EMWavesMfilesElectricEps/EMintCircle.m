clear all; clc; close all; clear workspace
input = initEM();           
      a = input.a;        
   gam0 = input.gamma_0; 
eps_sct = input.eps_sct;              
 mu_sct = input.mu_sct;
gam_sct = gam0 * sqrt(eps_sct*mu_sct); 
 Z_sct  = sqrt(mu_sct/eps_sct);  

% (1) Transform Cartesian coordinates to polar c0ordinates ----------------  
  xS = input.xS;   
  X1 = input.X1;  X2 = input.X2;
  rS = sqrt(xS(1)^2+xS(2)^2);    phiS = atan2(xS(2),xS(1)); 
   R = sqrt(X1.^2+X2.^2+1e-16);  PHI  = atan2(X2,X1); 

% (2) Compute coefficients of Bessel series expansion ---------------------
  arg0 = gam0 * input.a;  
  args = gam_sct*input.a; 
  M = 20;                                   % increase M for more accuracy
  A = zeros(1,M); B = zeros(1,M);
  for m = 1 : M                     
    Ib0 = besseli(m,arg0);     dIb0 =  besseli(m+1,arg0) + m/arg0 * Ib0; 
    Ibs = besseli(m,args);     dIbs =  besseli(m+1,args) + m/args * Ibs; 
    Kb0 = besselk(m,arg0);     dKb0 = -besselk(m+1,arg0) + m/arg0 * Kb0;
    denominator = Z_sct * dIbs*Kb0 - dKb0*Ibs;
    A(m)  =     -(Z_sct * dIbs*Ib0 - dIb0*Ibs) / denominator;
    B(m)  = (Ib0 + A(m) * Kb0) / Ibs;  
  end
  
% (3) Compute interior wave field using M terms of Bessel series expansion
    Er = zeros(size(R));  Ephi = zeros(size(R));  ZH3 = zeros(size(R));
  for m = 1 : M 
    args = gam_sct*R;   Ibs = besseli(m,args); 
                       dIbs = besseli(m+1,args) + m./args .* Ibs; 
                        KbS = besselk(m,gam0*rS);             
    Er   = Er   + B(m)*2*m^2.* Ibs .* KbS .* cos(m*(PHI-phiS));
    Ephi = Ephi - B(m)*2*m  .*dIbs .* KbS .* sin(m*(PHI-phiS));          
  end % m_loop
  Er   = (      1/eps_sct) * 1/(2*pi) * Er./R ./rS;
  Ephi = (gam_sct/eps_sct) * 1/(2*pi) * Ephi  ./rS; 
  
% (4) Transform to Cartesian vectors --------------------------------------
  E1 = cos(PHI) .* Er  - sin(PHI) .* Ephi;
  E2 = sin(PHI) .* Er  + cos(PHI) .* Ephi;

% (5) Plot contrast sources
  w_E{1} = input.CHI_eps .* E1;
  w_E{2} = input.CHI_eps .* E2;
  
  plotContrastSourcewE(w_E,input);                save CONTRASTSOURCE w_E; 