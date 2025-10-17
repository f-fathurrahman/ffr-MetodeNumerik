clear all; clc; close all;
input = init();   
c_0   = input.c_0;  c_sct = input.c_sct;

gam0    = input.gamma_0;               
gam_sct = input.gamma_0 * c_0/c_sct;   

% Transform Cartesian coordinates to polar ccordinates  
  xS = input.xS; rS = sqrt(xS(1)^2+xS(2)^2);    phiS = atan2(xS(2),xS(1)); 
  X1 = input.X1; X2 = input.X2; R = sqrt(X1.^2+X2.^2); PHI = atan2(X2,X1);    
  
% (1) Compute incident wave in closed form --------------------------------
      DIS    = sqrt(rS^2+R.^2-2*rS*R.*cos(phiS-PHI));
      u_inc = 1/(2*pi) .* besselk(0,gam0*DIS);

% (2) Compute coefficients of series expansion ----------------------------
      arg0 = gam0 * input.a;  args = gam_sct*input.a; 
      M = 50;                               % increase M for more accuracy
      A = zeros(1,M+1); B = zeros(1,M+1);
      for m = 0 : M                     
      Ib0 = besseli(m,arg0);     dIb0 =  besseli(m+1,arg0) + m/arg0 * Ib0; 
      Ibs = besseli(m,args);     dIbs =  besseli(m+1,args) + m/args * Ibs; 
      Kb0 = besselk(m,arg0);     dKb0 = -besselk(m+1,arg0) + m/arg0 * Kb0;
      A(m+1) = - (gam_sct * dIbs*Ib0 - gam0 * dIb0*Ibs) ...
                /(gam_sct * dIbs*Kb0 - gam0 * dKb0*Ibs);
      B(m+1) = (Ib0 + A(m+1) * Kb0) / Ibs;  
  end

% (3) Compute exterior, interior and total fields -------------------------
      factor = 1/(2*pi) * besselk(0,gam0*rS); 
      u_rfl = A(1) * factor * besselk(0,gam0*R);  
      u_int = B(1) * factor * besseli(0,gam_sct*R);
      for m = 1 : M
        factor = (1/pi) * besselk(m,gam0*rS) .* cos(m*(phiS-PHI));
        u_rfl  = u_rfl + A(m+1) * factor .* besselk(m,gam0*R);
        u_int  = u_int + B(m+1) * factor .* besseli(m,gam_sct*R);
      end
      u_rfl(ceil(input.N1/2),ceil(input.N2/2))=0;% Avoid NaN number at R=0
      u = (u_rfl+u_inc) .* (R >= input.a) +  u_int .* (R < input.a);  
  
      plotWavefield(u_inc,u,input);